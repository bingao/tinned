use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{
    AdjointMap, AdjointMode, MatrixAdd, MatrixMul, Number, TimeEvolution, ZeroOperator,
};
use crate::internal::{intern_expr, join_mapped};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{
    NumberTolerance, downcast_from_arc, expression_error, generic_expression_error, is_expr_type,
    is_zero_expr,
};

// Exponential adjoint map (or conjugation operation in Lie algebra):
// exp(ad_{X})(Y) = exp(X)*Y*exp(-X) (`left_action` is `true`), or
// exp(ad_{-X})(Y) = exp(-X)*Y*exp(X) (`left_action` is `false`).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct ExpAdjointMap {
    generator: Arc<dyn Expr>,
    // Whether generator and its derivatives commute
    generator_derivative_commute: bool,
    // `target` is mostly used for zeorth order
    target: Arc<dyn Expr>,
    // Indicates whether `target` is built as the time-evolution `generator`
    is_time_evolution: bool,
    left_action: bool,
    max_commutator_order: u32,
    // BCH expansion of the (differentiated) exponential adjoint map up to the
    // order of differentiation, i.e. the order of `derivative`. The key of
    // `BTreeMap` is the order of commutators and the value contains
    // corresponding BCH expansion terms.
    bch_expansion: BTreeMap<u32, Vec<Arc<dyn Expr>>>,
    derivative: PertMultichain,
}

impl ExpAdjointMap {
    // `is_rotation` indicates whether we make `generator` as a rotation operator
    #[inline]
    fn make_generator(generator: Arc<dyn Expr>, is_rotation: Option<bool>) -> Arc<dyn Expr> {
        let is_rotation = is_rotation.unwrap_or(false);

        if is_rotation {
            MatrixMul::new(vec![Number::imaginary_unit(), generator])
                .expect("ExpAdjointMap::make_generator() failed to build rotation operator")
        } else {
            generator
        }
    }

    #[inline]
    pub fn builder(
        generator: Arc<dyn Expr>,
        target: Arc<dyn Expr>,
        generator_derivative_commute: Option<bool>,
        is_rotation: Option<bool>,
    ) -> ExpAdjointMapBuilder {
        ExpAdjointMapBuilder {
            generator: Self::make_generator(generator, is_rotation),
            generator_derivative_commute,
            target,
            is_time_evolution: false,
            left_action: None,
            max_commutator_order: None,
            bch_expansion: None,
            derivative: None,
        }
    }

    #[inline]
    pub fn builder_time_evolution(
        generator: Arc<dyn Expr>,
        is_forward: bool,
        generator_derivative_commute: Option<bool>,
        is_rotation: Option<bool>,
    ) -> ExpAdjointMapBuilder {
        let generator = Self::make_generator(generator, is_rotation);

        let target = match TimeEvolution::builder(generator.clone()).is_forward(is_forward).build()
        {
            Ok(e) => e,
            Err(e) => panic!("ExpAdjointMap::builder_time_evolution() encounters: {e}"),
        };

        ExpAdjointMapBuilder {
            generator,
            generator_derivative_commute,
            target,
            is_time_evolution: true,
            left_action: None,
            max_commutator_order: None,
            bch_expansion: None,
            derivative: None,
        }
    }

    // This function is only used inside this file
    #[inline]
    fn with_bch_expansion(
        &self,
        bch_expansion: BTreeMap<u32, Vec<Arc<dyn Expr>>>,
        derivative: Option<PertMultichain>,
    ) -> ExpAdjointMapBuilder {
        ExpAdjointMapBuilder {
            generator: self.generator.clone(),
            generator_derivative_commute: Some(self.generator_derivative_commute),
            target: self.target.clone(),
            is_time_evolution: self.is_time_evolution,
            left_action: Some(self.left_action),
            max_commutator_order: Some(self.max_commutator_order),
            bch_expansion: Some(bch_expansion),
            derivative,
        }
    }

    #[inline]
    pub fn generator(&self) -> &Arc<dyn Expr> {
        &self.generator
    }

    #[inline]
    pub fn generator_derivative_commute(&self) -> bool {
        self.generator_derivative_commute
    }

    #[inline]
    pub fn target(&self) -> &Arc<dyn Expr> {
        &self.target
    }

    #[inline]
    pub fn is_time_evolution(&self) -> bool {
        self.is_time_evolution
    }

    #[inline]
    pub fn left_action(&self) -> bool {
        self.left_action
    }

    #[inline]
    pub fn max_commutator_order(&self) -> u32 {
        self.max_commutator_order
    }

    #[inline]
    pub fn bch_expansion(&self) -> &BTreeMap<u32, Vec<Arc<dyn Expr>>> {
        &self.bch_expansion
    }

    #[inline]
    pub fn derivative(&self) -> &PertMultichain {
        &self.derivative
    }

    // Rebuilds BCH expansion by applying a fallible operation to BCH expansion terms
    fn try_map_bch_expansion(
        &self,
        operation: impl Fn(&Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
        message: &'static str,
    ) -> Result<BTreeMap<u32, Vec<Arc<dyn Expr>>>, TinnedError> {
        let mut new_bch_expansion: BTreeMap<u32, Vec<Arc<dyn Expr>>> = BTreeMap::new();

        for (&order, terms) in &self.bch_expansion {
            let mut new_terms = Vec::with_capacity(terms.len());

            for expr in terms {
                let new_expr = operation(expr).map_err(|e| {
                    generic_expression_error(
                        format!("{} for BCH expansion term {}", message, expr),
                        self,
                        Some(Box::new(e)),
                    )
                })?;

                if !is_expr_type::<ZeroOperator>(&new_expr) {
                    new_terms.push(new_expr);
                }
            }

            if !new_terms.is_empty() {
                new_terms.sort();
                new_bch_expansion.insert(order, new_terms);
            }
        }

        Ok(new_bch_expansion)
    }

    // Applies a fallible operation to generator and BCH expansion terms
    fn transform_generator_and_bch_expansion(
        &self,
        operation: impl Fn(&Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
        message: &'static str,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let new_generator = operation(&self.generator).map_err(|e| {
            generic_expression_error(
                format!("{} for generator {}", message, &self.generator),
                self,
                Some(Box::new(e)),
            )
        })?;

        if is_zero_expr(&new_generator, None) {
            return Ok(ZeroOperator::new());
        }

        if self.bch_expansion.is_empty() {
            return if &new_generator == &self.generator {
                Ok(self.clone_expr())
            } else {
                self.with_bch_expansion(self.bch_expansion.clone(), Some(self.derivative.clone()))
                    .generator(new_generator)
                    .build()
            };
        }

        let new_bch_expansion = self.try_map_bch_expansion(operation, message)?;

        if new_bch_expansion.is_empty() {
            Ok(ZeroOperator::new())
        } else if new_bch_expansion == self.bch_expansion && &new_generator == &self.generator {
            Ok(self.clone_expr())
        } else {
            self.with_bch_expansion(new_bch_expansion, Some(self.derivative.clone()))
                .generator(new_generator)
                .build()
        }
    }

    // Applies a fallible operation to BCH expansion terms
    fn transform_bch_expansion(
        &self,
        operation: impl Fn(&Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
        message: &'static str,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let new_bch_expansion = self.try_map_bch_expansion(operation, message)?;

        if new_bch_expansion.is_empty() {
            Ok(ZeroOperator::new())
        } else if new_bch_expansion == self.bch_expansion {
            Ok(self.clone_expr())
        } else {
            self.with_bch_expansion(new_bch_expansion, Some(self.derivative.clone())).build()
        }
    }
}

#[derive(Debug)]
pub struct ExpAdjointMapBuilder {
    generator: Arc<dyn Expr>,
    generator_derivative_commute: Option<bool>,
    target: Arc<dyn Expr>,
    is_time_evolution: bool,
    left_action: Option<bool>,
    max_commutator_order: Option<u32>,
    bch_expansion: Option<BTreeMap<u32, Vec<Arc<dyn Expr>>>>,
    derivative: Option<PertMultichain>,
}

impl ExpAdjointMapBuilder {
    #[inline]
    fn generator(mut self, generator: Arc<dyn Expr>) -> Self {
        self.generator = generator;
        self
    }

    #[inline]
    pub fn generator_derivative_commute(mut self, generator_derivative_commute: bool) -> Self {
        self.generator_derivative_commute = Some(generator_derivative_commute);
        self
    }

    //#[inline]
    //fn target(mut self, target: Arc<dyn Expr>) -> Self {
    //    self.target = target;
    //    self
    //}

    #[inline]
    pub fn left_action(mut self, left_action: bool) -> Self {
        self.left_action = Some(left_action);
        self
    }

    #[inline]
    pub fn max_commutator_order(mut self, max_commutator_order: u32) -> Self {
        self.max_commutator_order = Some(max_commutator_order);
        self
    }

    //#[inline]
    //fn bch_expansion(mut self, bch_expansion: BTreeMap<u32, Vec<Arc<dyn Expr>>>) -> Self {
    //    self.bch_expansion = Some(bch_expansion);
    //    self
    //}

    //#[inline]
    //fn derivative(mut self, derivative: PertMultichain) -> Self {
    //    self.derivative = Some(derivative);
    //    self
    //}

    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.generator.is_scalar() {
            return Err(expression_error(
                "ExpAdjointMapBuilder::build() gets a scalar generator",
                &self.generator,
                None,
            ));
        }

        if self.target.is_scalar() {
            return Err(expression_error(
                "ExpAdjointMapBuilder::build() gets a scalar target",
                &self.target,
                None,
            ));
        }

        if is_expr_type::<ZeroOperator>(&self.generator)
            || is_expr_type::<ZeroOperator>(&self.target)
        {
            return Ok(self.target);
        }

        let generator_derivative_commute = self.generator_derivative_commute.unwrap_or(true);
        let left_action = self.left_action.unwrap_or(true);

        let max_commutator_order = self.max_commutator_order.unwrap_or(u32::MAX);

        if self.generator.has_unperturbed_term() && max_commutator_order == u32::MAX {
            return Err(expression_error(
                "ExpAdjointMapBuilder::build() gets a non-perturbing generator with infinite commutator order",
                &self.generator,
                None,
            ));
        }

        // Undifferentiated expression of exponential adjoint map is simply `target`
        let bch_expansion = self.bch_expansion.unwrap_or(if self.is_time_evolution {
            BTreeMap::new()
        } else {
            BTreeMap::from([(0, vec![self.target.clone()])])
        });

        let derivative = self.derivative.unwrap_or(PertMultichain::new());

        Ok(intern_expr(Arc::new(ExpAdjointMap {
            generator: self.generator,
            generator_derivative_commute,
            target: self.target,
            is_time_evolution: self.is_time_evolution,
            left_action,
            max_commutator_order,
            bch_expansion,
            derivative,
        })))
    }
}

impl ExprInternal for ExpAdjointMap {
    impl_expr_internal_methods!(ExpAdjointMap, true);

    #[inline]
    fn replace_one_in_children(
        &self,
        expr: &Arc<dyn Expr>,
        replacement: Arc<dyn Expr>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.bch_expansion.is_empty() {
            return Ok(self.clone_expr());
        }

        self.transform_bch_expansion(
            |term: &Arc<dyn Expr>| term.replace_one(expr, replacement.clone(), include_derivatives),
            "ExpAdjointMap::replace_one_in_children() failed",
        )
    }

    #[inline]
    fn replace_all_in_children(
        &self,
        map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.bch_expansion.is_empty() {
            return Ok(self.clone_expr());
        }

        self.transform_bch_expansion(
            |term: &Arc<dyn Expr>| term.replace_all(map, include_derivatives),
            "ExpAdjointMap::replace_all_in_children() failed",
        )
    }

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "ExpAdjointMap({}; {}; {}; {}; {}; {}; {{{}}}; [{}])",
            self.left_action,
            self.max_commutator_order,
            self.generator.hash_key(),
            self.generator_derivative_commute,
            self.target.hash_key(),
            self.is_time_evolution,
            join_mapped(self.bch_expansion.iter(), ";", |&(order, terms)| {
                format!("{}: [{}]", order, join_mapped(terms, ";", |term| (*term).hash_key()))
            }),
            self.derivative.hash_key(),
        )
    }

    #[inline]
    fn expr_order(&self) -> u32 {
        self.derivative.total_order()
    }

    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<ExpAdjointMap>(other) {
            // We treat exponential adjoint maps with different `left_action`'s equally
            self.max_commutator_order == op.max_commutator_order
                && self.generator.deep_eq_superchains(&op.generator)
                && self.generator_derivative_commute == op.generator_derivative_commute
                && self.target.deep_eq_superchains(&op.target)
                && self.is_time_evolution == op.is_time_evolution
                && self.derivative.is_subchain(&op.derivative)
        } else {
            false
        }
    }

    #[inline]
    fn eq_by_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<ExpAdjointMap>(other) {
            self.left_action == op.left_action
                && self.max_commutator_order == op.max_commutator_order
                && &self.generator == &op.generator
                && self.generator_derivative_commute == op.generator_derivative_commute
                && &self.target == &op.target
                && self.is_time_evolution == op.is_time_evolution
                && self.derivative.is_subchain(&op.derivative)
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for ExpAdjointMap {
    impl_expr_common_methods!(false);

    #[inline]
    fn has_unperturbed_term(&self) -> bool {
        // See equation (35), J. Comput. Chem. 45, 2136-2152 (2024).
        self.target.has_unperturbed_term()
    }

    //FIXME: test this function
    //FIXME: add function to do BCH expansion for self.max_commutator_order < u32::MAX
    fn substitute_zero_perturbations(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.bch_expansion.is_empty() {
            return Ok(ZeroOperator::new());
        }

        let new_bch_expansion = self.try_map_bch_expansion(
            |term| term.substitute_zero_perturbations(freq_tol.clone()),
            "ExpAdjointMap::substitute_zero_perturbations() failed",
        )?;

        if new_bch_expansion.is_empty() {
            return Ok(ZeroOperator::new());
        }

        // For `generator` as a perturbing operator, the
        // Baker-Campbell-Hausdorff (BCH) expansion is simply the BCH expansion
        // at zero perturbation strength.
        if !self.generator.has_unperturbed_term() {
            let terms = new_bch_expansion.into_values().flatten().collect();
            return MatrixAdd::new(terms);
        }

        // We need to apply `substitute_zero_perturbations()` for `generator`
        // and use its result for the exponential adjoint map at zero
        // perturbation strength
        let new_generator =
            self.generator.substitute_zero_perturbations(freq_tol).map_err(|e| {
                generic_expression_error(
                    "ExpAdjointMap::substitute_zero_perturbations() failed for generator",
                    self,
                    Some(Box::new(e)),
                )
            })?;

        self.with_bch_expansion(new_bch_expansion, Some(self.derivative.clone()))
            .generator(new_generator.clone())
            .build()
    }

    //FIXME: test this function
    // Equation (XX), ...
    fn differentiate(&self, s: Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        // For `target` as a time-evolution operator +/-i*d/dt, the first-order
        // derivative must be performed on the operator.
        if self.is_time_evolution && self.derivative.is_empty() {
            let diff_target = self.target.differentiate(s.clone()).map_err(|e| {
                generic_expression_error(
                    "ExpAdjointMap::differentiate() failed for target",
                    self,
                    Some(Box::new(e)),
                )
            })?;

            return if is_expr_type::<ZeroOperator>(&diff_target) {
                Ok(diff_target)
            } else {
                let slice: &[Arc<Perturbation>] = std::slice::from_ref(&s);
                self.with_bch_expansion(
                    BTreeMap::from([(1, vec![diff_target])]),
                    Some(PertMultichain::from_slice(slice)),
                )
                .build()
            };
        }

        // The next higher-order derivative of the exponential adjoint map is
        // computed from two ways. First, we differentiate each BCH expansion
        // term with respect to the perturbation `s`.
        let mut diff_bch_expansion = self.try_map_bch_expansion(
            |term| term.differentiate(s.clone()),
            "ExpAdjointMap::differentiate() failed",
        )?;

        // Secondly, for each BCH expansion term of the lower-order derivative,
        // we can build a commutator with generator as the differentiated
        // `generator `with respect to the perturbation `s`, and target as the
        // BCH expansion term.
        let diff_generator = self.generator.differentiate(s.clone()).map_err(|e| {
            generic_expression_error(
                "ExpAdjointMap::differentiate() failed for generator",
                self,
                Some(Box::new(e)),
            )
        })?;

        let new_deriv = Some(self.derivative.with_added_perturbation(s));

        if is_expr_type::<ZeroOperator>(&diff_generator) {
            return if diff_bch_expansion.is_empty() {
                Ok(ZeroOperator::new())
            } else {
                self.with_bch_expansion(diff_bch_expansion, new_deriv).build()
            };
        }

        // `adj_maps` contains commutators whose orders greater than the
        // maximum commutator order and should be extracted from the
        // exponential adjoint map
        let mut adj_maps = Vec::new();

        let adjoint_mode = if self.generator_derivative_commute {
            Some(AdjointMode::Commutative)
        } else {
            Some(AdjointMode::Symmetrized)
        };

        for (&order, terms) in &self.bch_expansion {
            if order < self.max_commutator_order {
                let diff_terms = diff_bch_expansion.entry(order + 1).or_default();
                diff_terms.reserve(terms.len());

                for expr in terms {
                    diff_terms.push(AdjointMap::new(
                        vec![diff_generator.clone()],
                        expr.clone(),
                        Some(self.left_action),
                        adjoint_mode,
                    )?);
                }

                diff_terms.sort();
            } else {
                adj_maps.reserve(terms.len());
                adj_maps.extend(terms.iter().cloned());
            }
        }

        let new_ead_map = self.with_bch_expansion(diff_bch_expansion, new_deriv).build()?;

        if adj_maps.is_empty() {
            Ok(new_ead_map)
        } else {
            adj_maps.push(new_ead_map);
            MatrixAdd::new(adj_maps)
        }
    }

    #[inline]
    fn eliminate(
        &self,
        parameter: Arc<dyn Expr>,
        perturbations: &[Arc<Perturbation>],
        min_order: u32,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        self.transform_generator_and_bch_expansion(
            |arg: &Arc<dyn Expr>| arg.eliminate(parameter.clone(), perturbations, min_order),
            "ExpAdjointMap::eliminate() failed",
        )
    }

    #[inline]
    fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
        if self.deep_eq_superchains(s) {
            BTreeMap::from([(self.expr_order(), HashSet::from([self.clone_expr()]))])
        } else {
            let mut result = self.generator.find_all(s);

            for terms in self.bch_expansion.values() {
                for term in terms {
                    for (order, subset) in term.find_all(s) {
                        result.entry(order).or_default().extend(subset);
                    }
                }
            }

            result
        }
    }

    #[inline]
    fn match_one(&self, s: &Arc<dyn Expr>, include_derivatives: bool) -> bool {
        // `target` should exists in `bch_expansion`
        self.match_one_self(s, include_derivatives)
            || self.generator.match_one(s, include_derivatives)
            || self
                .bch_expansion
                .values()
                .any(|terms| terms.iter().any(|expr| expr.match_one(s, include_derivatives)))
    }

    #[inline]
    fn match_any(&self, set: &HashSet<Arc<dyn Expr>>, include_derivatives: bool) -> bool {
        // `target` should exists in `bch_expansion`
        self.match_any_self(set, include_derivatives)
            || self.generator.match_any(set, include_derivatives)
            || self
                .bch_expansion
                .values()
                .any(|terms| terms.iter().any(|expr| expr.match_any(set, include_derivatives)))
    }

    #[inline]
    fn remove_one(&self, s: &Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_one_self(s, false) {
            return Ok(ZeroOperator::new());
        }

        self.transform_generator_and_bch_expansion(
            |arg: &Arc<dyn Expr>| arg.remove_one(s),
            "ExpAdjointMap::remove_one() failed",
        )
    }

    #[inline]
    fn remove_all(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_any_self(set, false) {
            return Ok(ZeroOperator::new());
        }

        self.transform_generator_and_bch_expansion(
            |arg: &Arc<dyn Expr>| arg.remove_all(set),
            "ExpAdjointMap::remove_all() failed",
        )
    }

    #[inline]
    fn retain_one(
        &self,
        s: &Arc<dyn Expr>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_one_self(s, include_derivatives) {
            return Ok(self.clone_expr());
        }

        if self.bch_expansion.is_empty() {
            return Ok(ZeroOperator::new());
        }

        self.transform_bch_expansion(
            |term: &Arc<dyn Expr>| term.retain_one(s, include_derivatives),
            "ExpAdjointMap::retain_one() failed",
        )
    }

    #[inline]
    fn retain_any(
        &self,
        set: &HashSet<Arc<dyn Expr>>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_any_self(set, include_derivatives) {
            return Ok(self.clone_expr());
        }

        if self.bch_expansion.is_empty() {
            return Ok(ZeroOperator::new());
        }

        self.transform_bch_expansion(
            |term: &Arc<dyn Expr>| term.retain_any(set, include_derivatives),
            "ExpAdjointMap::retain_any() failed",
        )
    }
}

impl PartialEq for ExpAdjointMap {
    fn eq(&self, other: &Self) -> bool {
        // We also compare `bch_expansion`, which may change after `substitute_zero_perturbations()`
        &self.generator == &other.generator
            && self.generator_derivative_commute == other.generator_derivative_commute
            && &self.target == &other.target
            && self.is_time_evolution == other.is_time_evolution
            && self.left_action == other.left_action
            && self.max_commutator_order == other.max_commutator_order
            && self.derivative == other.derivative
            && self.bch_expansion == other.bch_expansion
    }
}

impl Eq for ExpAdjointMap {}

impl std::fmt::Display for ExpAdjointMap {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.left_action {
            write!(f, "exp(ad[{}; {}", self.generator, self.generator_derivative_commute)?;
        } else {
            write!(f, "exp(ad[-({}); {}", self.generator, self.generator_derivative_commute)?;
        }

        if self.max_commutator_order < u32::MAX {
            write!(f, "; {}", self.max_commutator_order)?;
        }

        write!(f, "])({})^{}", self.target, self.derivative)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    test_struct_safety!(ExpAdjointMap);
}
