use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{MatrixAdd, ZeroOperator};
use crate::internal::{intern_expr, join_mapped, sort_expressions_grouped_by};
use crate::perturbations::Perturbation;
use crate::public::{
    NumberTolerance, downcast_from_arc, expression_error, generic_expression_error, is_expr_type,
    is_zero_expr,
};

// Mode of an adjoint map or its generators
#[cfg_attr(feature = "ffi", safer_ffi::derive_ReprC)]
#[repr(u32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum AdjointMode {
    Commutative, // generators and their derivatives commute -> canonical order
    Symmetrized, // symmetrized nested commutator
    Ordered,     // noncommutative generators, order preserved
}

impl std::fmt::Display for AdjointMode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let s = match self {
            AdjointMode::Commutative => "commutative",
            AdjointMode::Symmetrized => "symmetrized",
            AdjointMode::Ordered => "ordered",
        };
        write!(f, "{s}")
    }
}

// Either adjoint map (or adjoint action, adjoint representation)
// [xn, [..., [x1, [x0, y]]...]], or its "right" acting version
// [[...[[y, x0], x1], ...], xn]. Users decide which version they are using.
// The default is acting from the left side.
//
// `generators` holds x0, x1, ..., xn, and y is stored in `target`.
//
// For `AdjointMode` as `Commutative` and `Symmetrized`, we sort `generators`
// in a deterministic way without changing the result of adjoint map.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct AdjointMap {
    generators: Vec<Arc<dyn Expr>>,
    target: Arc<dyn Expr>,
    left_action: bool,
    adjoint_mode: AdjointMode,
}

impl AdjointMap {
    pub fn new(
        generators: Vec<Arc<dyn Expr>>,
        target: Arc<dyn Expr>,
        left_action: Option<bool>,
        adjoint_mode: Option<AdjointMode>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if target.is_scalar() {
            return Err(expression_error("AdjointMap::new() gets a scalar target", &target, None));
        }
        if is_expr_type::<ZeroOperator>(&target) {
            return Ok(target);
        }

        for generator in &generators {
            if generator.is_scalar() {
                return Err(expression_error(
                    "AdjointMap::new() got a scalar generator",
                    generator,
                    None,
                ));
            }
            if is_expr_type::<ZeroOperator>(&generator) {
                return Ok(ZeroOperator::new());
            }
        }

        let left_action = left_action.unwrap_or(true);
        let adjoint_mode = adjoint_mode.unwrap_or(AdjointMode::Commutative);

        let generators = Self::canonicalize_generators(generators, adjoint_mode);

        Ok(intern_expr(Arc::new(Self {
            generators,
            target,
            left_action,
            adjoint_mode,
        })))
    }

    #[inline]
    fn canonicalize_generators(
        generators: Vec<Arc<dyn Expr>>,
        adjoint_mode: AdjointMode,
    ) -> Vec<Arc<dyn Expr>> {
        // Sort generators according to `expr_order()`
        if matches!(adjoint_mode, AdjointMode::Ordered) {
            generators
        } else {
            sort_expressions_grouped_by(&generators, |e| e.expr_order())
        }
    }

    #[inline]
    fn with_new_generators(
        &self,
        generators: Vec<Arc<dyn Expr>>,
        adjoint_mode: Option<AdjointMode>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let adjoint_mode = adjoint_mode.unwrap_or(self.adjoint_mode);

        Self::new(generators, self.target.clone(), Some(self.left_action), Some(adjoint_mode))
    }

    #[inline]
    fn with_new_target(&self, target: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if target.is_scalar() {
            return Err(expression_error(
                "AdjointMap::with_new_target() gets a scalar target",
                &target,
                None,
            ));
        }
        if is_expr_type::<ZeroOperator>(&target) {
            return Ok(target);
        }

        Ok(intern_expr(Arc::new(Self {
            generators: self.generators.clone(),
            target,
            left_action: self.left_action,
            adjoint_mode: self.adjoint_mode,
        })))
    }

    #[inline]
    pub(crate) fn with_added_generator(
        &self,
        generator: Arc<dyn Expr>,
        adjoint_mode: Option<AdjointMode>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let mut generators = self.generators.clone();
        generators.push(generator);

        let adjoint_mode = adjoint_mode.unwrap_or(self.adjoint_mode);

        Self::new(generators, self.target.clone(), Some(self.left_action), Some(adjoint_mode))
    }

    #[inline]
    pub fn generators(&self) -> &[Arc<dyn Expr>] {
        &self.generators
    }

    #[inline]
    pub fn target(&self) -> &Arc<dyn Expr> {
        &self.target
    }

    #[inline]
    pub fn left_action(&self) -> bool {
        self.left_action
    }

    #[inline]
    pub fn adjoint_mode(&self) -> AdjointMode {
        self.adjoint_mode
    }

    // Checks modes of two adjoint maps, and returns comparable generators if
    // modes are compatible
    fn comparable_generators_if_compatible(
        &self,
        other: &Self,
    ) -> Option<(Vec<Arc<dyn Expr>>, Vec<Arc<dyn Expr>>)> {
        let same_mode = self.adjoint_mode == other.adjoint_mode;
        let either_symmetrized = matches!(self.adjoint_mode, AdjointMode::Symmetrized)
            || matches!(other.adjoint_mode, AdjointMode::Symmetrized);

        // An adjoint map with `Symmetrized` mode can only be compared to the
        // other with the same mode
        if either_symmetrized && !same_mode {
            return None;
        }

        // If one adjoint map has mode `Commutative`, it may be equal to the
        // other even with `Ordered` mode
        let self_generators = if matches!(self.adjoint_mode, AdjointMode::Ordered) && !same_mode {
            Self::canonicalize_generators(self.generators.clone(), AdjointMode::Commutative)
        } else {
            self.generators.clone()
        };

        let other_generators = if matches!(other.adjoint_mode, AdjointMode::Ordered) && !same_mode {
            Self::canonicalize_generators(other.generators.clone(), AdjointMode::Commutative)
        } else {
            other.generators.clone()
        };

        Some((self_generators, other_generators))
    }

    // Applies a fallible operation to the target and generators of this
    // AdjointMap and reconstructs the expression if needed.
    //
    // Semantics:
    // - The provided operation is applied to the target and each generator.
    // - If the resulting target is zero, the entire AdjointMap collapses to ZeroOperator.
    // - If any resulting generator is zero, the entire AdjointMap collapses to ZeroOperator.
    // - If no child changes, the original expression is returned.
    // - Otherwise, a new AdjointMap is constructed with the transformed children.
    fn transform_children_any_zero(
        &self,
        operation: impl Fn(&Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
        message: &'static str,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let new_target = operation(&self.target).map_err(|e| {
            generic_expression_error(format!("{} for target", message), self, Some(Box::new(e)))
        })?;

        if is_zero_expr(&new_target, None) {
            return Ok(ZeroOperator::new());
        }

        let mut changed = !Arc::ptr_eq(&new_target, &self.target) && &new_target != &self.target;

        let mut new_generators = Vec::with_capacity(self.generators.len());

        for generator in &self.generators {
            let new_generator = operation(generator).map_err(|e| {
                generic_expression_error(
                    format!("{} for generator {}", message, generator),
                    self,
                    Some(Box::new(e)),
                )
            })?;

            if is_zero_expr(&new_generator, None) {
                return Ok(ZeroOperator::new());
            }

            if !changed && !Arc::ptr_eq(&new_generator, generator) && &new_generator != generator {
                changed = true;
            }

            new_generators.push(new_generator);
        }

        if changed {
            Self::new(new_generators, new_target, Some(self.left_action), Some(self.adjoint_mode))
        } else {
            Ok(self.clone_expr())
        }
    }

    // Applies a fallible operation to the target and generators of this
    // AdjointMap and reconstructs the expression if needed.
    //
    // Semantics:
    // - The provided operation is applied to the target and each generator.
    // - If the resulting target and all resulting generators are zero, the entire AdjointMap collapses to ZeroOperator.
    // - If only some resulting children are zero, those children are treated as
    //   unchanged and their original child expressions are retained.
    // - If no child changes, the original expression is returned.
    // - Otherwise, a new AdjointMap is constructed with the transformed children.
    fn transform_children_all_zero(
        &self,
        match_self: impl FnOnce(&Self, bool) -> bool,
        include_derivatives: bool,
        retain_child: impl Fn(&Arc<dyn Expr>, bool) -> Result<Arc<dyn Expr>, TinnedError>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if match_self(self, include_derivatives) {
            return Ok(self.clone_expr());
        }

        let retained_target = retain_child(&self.target, include_derivatives).map_err(|e| {
            generic_expression_error(
                "AdjointMap::transform_children_all_zero() failed for target",
                self,
                Some(Box::new(e)),
            )
        })?;

        let target_is_zero = is_zero_expr(&retained_target, None);
        let new_target = if target_is_zero {
            self.target.clone()
        } else {
            retained_target
        };

        let mut changed = !Arc::ptr_eq(&new_target, &self.target) && &new_target != &self.target;

        let mut all_zero = target_is_zero;
        let mut new_generators = Vec::with_capacity(self.generators.len());

        for generator in &self.generators {
            let retained_generator = retain_child(generator, include_derivatives).map_err(|e| {
                generic_expression_error(
                    format!(
                        "AdjointMap::transform_children_all_zero() failed for generator {}",
                        generator
                    ),
                    self,
                    Some(Box::new(e)),
                )
            })?;

            if is_zero_expr(&retained_generator, None) {
                new_generators.push(generator.clone());
            } else {
                all_zero = false;

                if !changed
                    && !Arc::ptr_eq(&retained_generator, generator)
                    && &retained_generator != generator
                {
                    changed = true;
                }

                new_generators.push(retained_generator);
            }
        }

        if all_zero {
            Ok(ZeroOperator::new())
        } else if changed {
            Self::new(new_generators, new_target, Some(self.left_action), Some(self.adjoint_mode))
        } else {
            Ok(self.clone_expr())
        }
    }
}

impl ExprInternal for AdjointMap {
    impl_expr_internal_methods!(AdjointMap, false);

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "AdjointMap([{}]; {}; {}; {})",
            join_mapped(&self.generators, ";", |generator| generator.hash_key()),
            self.target.hash_key(),
            self.left_action,
            self.adjoint_mode,
        )
    }

    #[inline]
    fn expr_order(&self) -> u32 {
        let mut expr_order = self.target.expr_order();

        for generator in &self.generators {
            expr_order += generator.expr_order()
        }

        expr_order
    }

    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        let Some(op) = downcast_from_arc::<AdjointMap>(other) else {
            return false;
        };

        // We treat adjoint maps with different `left_action`'s equally
        if !self.target.deep_eq_superchains(&op.target)
            || self.generators.len() != op.generators.len()
        {
            return false;
        }

        let Some((self_generators, other_generators)) =
            self.comparable_generators_if_compatible(&op)
        else {
            return false;
        };

        self_generators.iter().zip(&other_generators).all(|(a, b)| a.deep_eq_superchains(b))
    }

    fn replace_one_in_children(
        &self,
        expr: &Arc<dyn Expr>,
        replacement: Arc<dyn Expr>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        self.transform_children_any_zero(
            |x| x.replace_one(expr, replacement.clone(), include_derivatives),
            "AdjointMap::replace_one_in_children() failed",
        )
    }

    fn replace_all_in_children(
        &self,
        map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        self.transform_children_any_zero(
            |x| x.replace_all(map, include_derivatives),
            "AdjointMap::replace_all_in_children() failed",
        )
    }
}

#[typetag::serde]
impl Expr for AdjointMap {
    impl_expr_common_methods!(false);

    #[inline]
    fn has_unperturbed_term(&self) -> bool {
        self.target.has_unperturbed_term()
            && self.generators.iter().all(|g| g.has_unperturbed_term())
    }

    fn substitute_zero_perturbations(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        self.transform_children_any_zero(
            |x| x.substitute_zero_perturbations(freq_tol.clone()),
            "AdjointMap::substitute_zero_perturbations() failed",
        )
    }

    fn differentiate(&self, s: Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        // Precompute the derivative of each x and store it
        let with_context = |f: &Arc<dyn Expr>| {
            f.differentiate(s.clone()).map_err(|e| {
                generic_expression_error(
                    "AdjointMap::differentiate() failed for generators",
                    self,
                    Some(Box::new(e)),
                )
            })
        };

        let diff_generators: Vec<Arc<dyn Expr>> =
            self.generators.iter().map(with_context).collect::<Result<_, _>>()?;

        let mut results = Vec::with_capacity(diff_generators.len() + 1);

        for (i, diff) in diff_generators.iter().enumerate() {
            // Skip derivative = 0 to avoid 0 * [...] = 0
            if is_zero_expr(diff, None) {
                continue;
            }

            let mut new_generators = self.generators.clone();
            // For each x, replace it with its derivative while keeping others
            // intact
            new_generators[i] = diff.clone();

            results.push(self.with_new_generators(new_generators, Some(self.adjoint_mode))?);
        }

        let diff_target = self.target.differentiate(s).map_err(|e| {
            generic_expression_error(
                "AdjointMap::differentiate() failed for target",
                self,
                Some(Box::new(e)),
            )
        })?;

        if !is_zero_expr(&diff_target, None) {
            results.push(self.with_new_target(diff_target)?);
        }

        MatrixAdd::new(results)
    }

    fn eliminate(
        &self,
        parameter: Arc<dyn Expr>,
        perturbations: &[Arc<Perturbation>],
        min_order: u32,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        self.transform_children_any_zero(
            |x| x.eliminate(parameter.clone(), perturbations, min_order),
            "AdjointMap::eliminate() failed",
        )
    }

    #[inline]
    fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
        if self.deep_eq_superchains(s) {
            return BTreeMap::from([(self.expr_order(), HashSet::from([self.clone_expr()]))]);
        }

        let mut result = self.target.find_all(s);

        for x in &self.generators {
            for (order, subset) in x.find_all(s) {
                result.entry(order).or_default().extend(subset);
            }
        }

        result
    }

    #[inline]
    fn match_one(&self, s: &Arc<dyn Expr>, include_derivatives: bool) -> bool {
        self.match_one_self(s, include_derivatives)
            || self.target.match_one(s, include_derivatives)
            || self.generators.iter().any(|x| x.match_one(s, include_derivatives))
    }

    #[inline]
    fn match_any(&self, set: &HashSet<Arc<dyn Expr>>, include_derivatives: bool) -> bool {
        self.match_any_self(set, include_derivatives)
            || self.target.match_any(set, include_derivatives)
            || self.generators.iter().any(|x| x.match_any(set, include_derivatives))
    }

    fn remove_one(&self, s: &Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_one_self(s, false) {
            return Ok(ZeroOperator::new());
        }

        self.transform_children_any_zero(|x| x.remove_one(s), "AdjointMap::remove_one() failed")
    }

    fn remove_all(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_any_self(set, false) {
            return Ok(ZeroOperator::new());
        }

        self.transform_children_any_zero(|x| x.remove_all(set), "AdjointMap::remove_all() failed")
    }

    fn retain_one(
        &self,
        s: &Arc<dyn Expr>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        self.transform_children_all_zero(
            |this, include_derivatives| this.match_one_self(s, include_derivatives),
            include_derivatives,
            |expr, include_derivatives| expr.retain_one(s, include_derivatives),
        )
    }

    fn retain_any(
        &self,
        set: &HashSet<Arc<dyn Expr>>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        self.transform_children_all_zero(
            |this, include_derivatives| this.match_any_self(set, include_derivatives),
            include_derivatives,
            |expr, include_derivatives| expr.retain_any(set, include_derivatives),
        )
    }
}

impl PartialEq for AdjointMap {
    fn eq(&self, other: &Self) -> bool {
        let len = self.generators.len();

        if &self.target != &other.target
            || (self.left_action != other.left_action && len % 2 == 1)
            || len != other.generators.len()
        {
            return false;
        }

        let Some((self_generators, other_generators)) =
            self.comparable_generators_if_compatible(other)
        else {
            return false;
        };

        self_generators == other_generators
    }
}

impl Eq for AdjointMap {}

impl std::fmt::Display for AdjointMap {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.left_action {
            write!(
                f,
                "[{},{}{}_{{{}}}",
                join_mapped(&self.generators, ",[", |generator| generator.to_string()),
                self.target,
                "]".repeat(self.generators.len()),
                self.adjoint_mode
            )
        } else {
            write!(
                f,
                "{}{},{}]_{{{}}}",
                "[".repeat(self.generators.len()),
                self.target,
                join_mapped(&self.generators, "],", |generator| generator.to_string()),
                self.adjoint_mode
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    test_struct_safety!(AdjointMap);
}
