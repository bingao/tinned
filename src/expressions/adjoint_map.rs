use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{MatrixAdd, ZeroOperator};
use crate::internal::{
    intern_expr, multi_expression_format, multi_expression_hash, sort_expressions_grouped_by,
};
use crate::perturbations::Perturbation;
use crate::public::{
    NumberTolerance, downcast_from_arc, downcast_from_ref, expression_error,
    generic_expression_error, is_expr_type, is_zero_expr,
};

// Either adjoint map (or adjoint action, adjoint representation)
// [xn, [..., [x1, [x0, y]]...]], or its "right" acting version
// [[...[[y, x0], x1], ...], xn]. Users decide which version they are using.
// The default is acting from the left side.
//
// `generators` holds x0, x1, ..., xn, and y is stored in `target`.
// All generators are commutative, i.e. [xi, xj] = 0 for all 0 <= i, j <= n.
// So, we can sort `generators` in some way without changing the result of
// adjoint map.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct AdjointMap {
    generators: Vec<Arc<dyn Expr>>,
    target: Arc<dyn Expr>,
    left_action: bool,
}

impl AdjointMap {
    pub fn new(
        generators: Vec<Arc<dyn Expr>>,
        target: Arc<dyn Expr>,
        left_action: Option<bool>,
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

        // Sort generators according to `total_order()`
        let sorted = sort_expressions_grouped_by(&generators, |e| e.total_order());

        let left_action = left_action.unwrap_or(true);

        Ok(intern_expr(Arc::new(Self {
            generators: sorted,
            target,
            left_action,
        })))
    }

    #[inline]
    fn with_new_generators(
        &self,
        generators: Vec<Arc<dyn Expr>>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Self::new(generators, self.target.clone(), Some(self.left_action))
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
        })))
    }

    #[inline]
    pub(crate) fn with_added_generator(
        &self,
        generator: Arc<dyn Expr>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let mut generators = self.generators.clone();
        generators.push(generator);

        Self::new(generators, self.target.clone(), Some(self.left_action))
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
}

impl ExprInternal for AdjointMap {
    impl_expr_internal_methods!(AdjointMap);

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "AdjointMap([{}]; {}; {})",
            multi_expression_hash(&self.generators, ";"),
            self.target.hash_key(),
            self.left_action,
        )
    }

    #[inline]
    fn match_for_find_all(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<AdjointMap>(other) {
            let len = self.generators.len();

            // We find adjoint maps with `left_action` either `true` or `false`
            if !self.target.match_for_find_all(&op.target) || len != op.generators.len() {
                return false;
            }

            self.generators.iter().zip(&op.generators).all(|(a, b)| a.match_for_find_all(b))
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for AdjointMap {
    impl_expr_common_methods!(false);

    #[inline]
    fn clean_temporum(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        impl_adjoint_map_operation!(
            self,
            |x: &Arc<dyn Expr>| x.clean_temporum(freq_tol.clone()),
            "AdjointMap::clean_temporum() failed"
        )
    }

    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        // Precompute the derivative of each x and store it
        let with_context = |f: &Arc<dyn Expr>| {
            f.differentiate(s).map_err(|e| {
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

            results.push(self.with_new_generators(new_generators)?);
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

    #[inline]
    fn eliminate(
        &self,
        parameter: &Arc<dyn Expr>,
        perturbations: &[Arc<Perturbation>],
        min_order: u32,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        impl_adjoint_map_operation!(
            self,
            |x: &Arc<dyn Expr>| x.eliminate(parameter, perturbations, min_order),
            "AdjointMap::eliminate() failed"
        )
    }

    #[inline]
    fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
        if self.generators.iter().any(|x| x.exist_any(set)) {
            return true;
        }

        set.iter().any(|expr| self.eq_expr(expr.as_ref())) || self.target.exist_any(set)
    }

    #[inline]
    fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
        if self.match_for_find_all(s) {
            return BTreeMap::from([(self.total_order(), HashSet::from([self.clone_expr()]))]);
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
    fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
            return Ok(ZeroOperator::new());
        }

        impl_adjoint_map_operation!(
            self,
            |x: &Arc<dyn Expr>| x.remove(set),
            "AdjointMap::remove() failed"
        )
    }

    #[inline]
    fn retain(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
            return Ok(self.clone_expr());
        }

        impl_adjoint_map_operation!(
            self,
            |x: &Arc<dyn Expr>| x.retain(set),
            "AdjointMap::retain() failed"
        )
    }

    #[inline]
    fn replace(
        &self,
        map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if let Some((_, value)) = map.iter().find(|(key, _)| self.eq_expr(key.as_ref())) {
            return Ok(value.clone());
        }

        impl_adjoint_map_operation!(
            self,
            |x: &Arc<dyn Expr>| x.replace(map),
            "AdjointMap::replace() failed"
        )
    }

    #[inline]
    fn replace_all(
        &self,
        map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if let Some((_, value)) = map.iter().find(|(key, _)| self.match_for_replace_all(key)) {
            return Ok(value.clone());
        }

        impl_adjoint_map_operation!(
            self,
            |x: &Arc<dyn Expr>| x.replace_all(map),
            "AdjointMap::replace_all() failed"
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

        self.generators == other.generators
    }
}

impl Eq for AdjointMap {}

impl std::fmt::Display for AdjointMap {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.left_action {
            write!(
                f,
                "[{},{}{}",
                multi_expression_format(&self.generators, ",["),
                self.target,
                "]".repeat(self.generators.len())
            )
        } else {
            write!(
                f,
                "{}{},{}]",
                "[".repeat(self.generators.len()),
                self.target,
                multi_expression_format(&self.generators, "],")
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    test_struct_safety!(AdjointMap);
}
