use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{
    AdjointMap, AdjointMode, MatrixAdd, MatrixMul, Number, TimeEvolution, ZeroOperator,
};
use crate::internal::intern_expr;
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{
    NumberTolerance, downcast_from_arc, expression_error, generic_expression_error, is_expr_type,
    is_zero_expr,
};
use crate::unreachable_error;

// Exponential adjoint map (or conjugation operation in Lie algebra):
// exp(ad_{X})(Y) = exp(X)*Y*exp(-X) (`left_action` is `true`), or
// exp(ad_{-X})(Y) = exp(-X)*Y*exp(X) (`left_action` is `false`).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct ExpAdjointMap {
    generator: Arc<dyn Expr>,
    // Whether generator and its derivatives commute
    generator_derivative_commute: bool,
    target: Arc<dyn Expr>,
    is_time_evolution: bool,
    left_action: bool,
    max_commutator_order: u32,
    // `result` contains differentiated expression of exponential adjoint map
    result: Arc<dyn Expr>,
    derivative: PertMultichain,
}

impl ExpAdjointMap {
    #[inline]
    pub fn builder(
        generator: Arc<dyn Expr>,
        target: Arc<dyn Expr>,
        generator_derivative_commute: Option<bool>,
    ) -> ExpAdjointMapBuilder {
        ExpAdjointMapBuilder {
            generator,
            generator_derivative_commute,
            target,
            is_time_evolution: false,
            left_action: None,
            max_commutator_order: None,
            result: None,
            derivative: None,
        }
    }

    #[inline]
    pub fn builder_time_evolution(
        generator: Arc<dyn Expr>,
        is_forward: bool,
        generator_derivative_commute: Option<bool>,
    ) -> ExpAdjointMapBuilder {
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
            result: None,
            derivative: None,
        }
    }

    // This function and `with_result_and_derivative()` are only used inside this file
    #[inline]
    fn with_result(&self, result: Arc<dyn Expr>) -> ExpAdjointMapBuilder {
        ExpAdjointMapBuilder {
            generator: self.generator.clone(),
            generator_derivative_commute: Some(self.generator_derivative_commute),
            target: self.target.clone(),
            is_time_evolution: self.is_time_evolution,
            left_action: Some(self.left_action),
            max_commutator_order: Some(self.max_commutator_order),
            result: Some(result),
            derivative: Some(self.derivative.clone()),
        }
    }

    #[inline]
    fn with_result_and_derivative(
        &self,
        result: Arc<dyn Expr>,
        derivative: PertMultichain,
    ) -> ExpAdjointMapBuilder {
        ExpAdjointMapBuilder {
            generator: self.generator.clone(),
            generator_derivative_commute: Some(self.generator_derivative_commute),
            target: self.target.clone(),
            is_time_evolution: self.is_time_evolution,
            left_action: Some(self.left_action),
            max_commutator_order: Some(self.max_commutator_order),
            result: Some(result),
            derivative: Some(derivative),
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
    pub fn result(&self) -> &Arc<dyn Expr> {
        &self.result
    }

    #[inline]
    pub fn derivative(&self) -> &PertMultichain {
        &self.derivative
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
    result: Option<Arc<dyn Expr>>,
    derivative: Option<PertMultichain>,
}

impl ExpAdjointMapBuilder {
    #[inline]
    pub fn generator_derivative_commute(mut self, generator_derivative_commute: bool) -> Self {
        self.generator_derivative_commute = Some(generator_derivative_commute);
        self
    }

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
    //fn result(mut self, result: Arc<dyn Expr>) -> Self {
    //    self.result = Some(result);
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
        let result = self.result.unwrap_or(self.target.clone());
        if is_expr_type::<ZeroOperator>(&result) {
            return Ok(ZeroOperator::new());
        }

        let derivative = self.derivative.unwrap_or(PertMultichain::new());

        Ok(intern_expr(Arc::new(ExpAdjointMap {
            generator: self.generator,
            generator_derivative_commute,
            target: self.target,
            is_time_evolution: self.is_time_evolution,
            left_action,
            max_commutator_order,
            result,
            derivative,
        })))
    }
}

impl ExprInternal for ExpAdjointMap {
    impl_unary_expr_internal_methods!(
        ExpAdjointMap,
        False,
        result,
        true,
        |this: &ExpAdjointMap, arg| { this.with_result(arg).build() }
    );

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "ExpAdjointMap({}; {}; {}; {}; {}; {}; {}; [{}])",
            self.left_action,
            self.max_commutator_order,
            self.generator.hash_key(),
            self.generator_derivative_commute,
            self.target.hash_key(),
            self.is_time_evolution,
            self.result.hash_key(),
            self.derivative.hash_key(),
        )
    }

    #[inline]
    fn total_order(&self) -> u32 {
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
    impl_unary_expr_common_methods!(ExpAdjointMap, False, result, |this: &ExpAdjointMap, arg| this
        .with_result(arg)
        .build());

    #[inline]
    fn has_unperturbed_term(&self) -> bool {
        // See equation (35), J. Comput. Chem. 45, 2136-2152 (2024).
        self.target.has_unperturbed_term()
    }

    //FIXME: test this function
    #[inline]
    fn substitute_zero_perturbations(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let result = self.result.substitute_zero_perturbations(freq_tol.clone()).map_err(|e| {
            generic_expression_error(
                "ExpAdjointMap::substitute_zero_perturbations() failed for result",
                self,
                Some(Box::new(e)),
            )
        })?;

        if is_zero_expr(&result, freq_tol.clone()) {
            return Ok(ZeroOperator::new());
        }

        // For finite commutator order like coupled-cluster theory, we return
        // `ExpAdjointMap` with updated `result`.
        if self.max_commutator_order < u32::MAX {
            return ExpAdjointMap::with_result(&self, result).build();
        }

        // For `generator` as a perturbing operator, the
        // Baker-Campbell-Hausdorff (BCH) expansion is simply `result` at zero
        // perturbation strength.
        if !self.generator.has_unperturbed_term() {
            return Ok(result);
        }

        // `ExpAdjointMapBuilder::build()` should prevent this error, but it is
        // worthy of checking again.
        if self.max_commutator_order == u32::MAX {
            return Err(unreachable_error(
                "ExpAdjointMap::substitute_zero_perturbations() gets a non-perturbing generator with infinite commutator order",
                &self.generator,
                None,
            ));
        }

        // We need to apply `substitute_zero_perturbations()` for `generator`
        // and use its result for the BCH expansion.
        let generator = self.generator.substitute_zero_perturbations(freq_tol).map_err(|e| {
            generic_expression_error(
                "ExpAdjointMap::substitute_zero_perturbations() failed for generator",
                self,
                Some(Box::new(e)),
            )
        })?;

        // This method expands the exponential adjoint map using the
        // Baker-Campbell-Hausdorff (BCH) expansion.
        #[inline]
        fn do_bch_expansion(
            generator: &Arc<dyn Expr>,
            target: &Arc<dyn Expr>,
            max_commutator_order: u32,
            left_action: Option<bool>,
            adjoint_mode: Option<AdjointMode>,
        ) -> Result<Vec<Arc<dyn Expr>>, TinnedError> {
            let mut terms = Vec::with_capacity((max_commutator_order as usize) + 1);
            terms.push(target.clone());

            if max_commutator_order == 0 {
                return Ok(terms);
            }

            let mut generators = Vec::with_capacity(max_commutator_order as usize);
            let mut denom: i64 = 1;

            for order in 1..=max_commutator_order {
                generators.push(generator.clone());

                let adj_map =
                    AdjointMap::new(generators.clone(), target.clone(), left_action, adjoint_mode)?;

                if order == 1 {
                    terms.push(adj_map);
                } else {
                    denom *= order as i64;
                    let coefficient =
                        Number::from_rational(num_rational::Rational64::new(1, denom));
                    terms.push(MatrixMul::new(vec![coefficient, adj_map])?);
                }
            }

            Ok(terms)
        }

        // Now, we will apply the BCH expansion for `generator` and `result`,
        // which will result into an `MatrixAdd`.  `result` is either (i)
        // undifferentiated `target`, (ii) an `AdjointMap` when only
        // `generator` was differentiated, or (iii) an `MatrixAdd` of
        // `AdjointMap`'s and differentiated `target`.
        let estimated_terms = if let Some(mat_add) = downcast_from_arc::<MatrixAdd>(&result) {
            mat_add.terms().len() * ((self.max_commutator_order as usize) + 1)
        } else {
            (self.max_commutator_order as usize) + 1
        };

        let mut adj_maps = Vec::with_capacity(estimated_terms);

        let adjoint_mode = if self.generator_derivative_commute {
            Some(AdjointMode::Commutative)
        } else {
            Some(AdjointMode::Symmetric)
        };

        // A helper closure to process `result` or its terms when it is an `MatrixAdd`
        let mut process_term = |term: &Arc<dyn Expr>| -> Result<(), TinnedError> {
            let commutator_order = if let Some(adj_map) = downcast_from_arc::<AdjointMap>(term) {
                adj_map.generators().len() as u32
            } else {
                0
            };

            if commutator_order > self.max_commutator_order {
                return Err(unreachable_error(
                    format!(
                        "ExpAdjointMap::substitute_zero_perturbations() got a target violating its maximum commutator order {}",
                        self.max_commutator_order
                    ),
                    term,
                    None,
                ));
            }

            let bch_terms = do_bch_expansion(
                &generator,
                term,
                self.max_commutator_order - commutator_order,
                Some(self.left_action),
                adjoint_mode,
            )?;
            adj_maps.extend(bch_terms);
            Ok(())
        };

        if let Some(mat_add) = downcast_from_arc::<MatrixAdd>(&result) {
            for term in mat_add.terms() {
                process_term(term)?;
            }
        } else {
            process_term(&result)?;
        }

        MatrixAdd::new(adj_maps)
    }

    //FIXME: test this function
    // Equation (47)
    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        // `result` is either (i) undifferentiated `target`, (ii) an
        // `AdjointMap` when only `generator` was differentiated, or (iii) an
        // `MatrixAdd` of `AdjointMap`'s and differentiated `target`. We first
        // differentiate `result` with respect to `s`, which gives us all
        // differentiated terms (including the differentiation on each
        // `AdjointMap`'s field `target`) with number of generators fixed for
        // each `AdjointMap`.
        let diff_result = self.result.differentiate(s).map_err(|e| {
            generic_expression_error(
                "ExpAdjointMap::differentiate() failed for result",
                self,
                Some(Box::new(e)),
            )
        })?;

        // The first order derivative must be performed on the time evolution
        // operator +/-i*d/dt
        if self.is_time_evolution && self.derivative.is_empty() {
            return if is_expr_type::<ZeroOperator>(&diff_result) {
                Ok(diff_result)
            } else {
                let slice: &[Arc<Perturbation>] = std::slice::from_ref(s);
                self.with_result_and_derivative(diff_result, PertMultichain::from_slice(slice))
                    .build()
            };
        }

        // For each previous differentiated `AdjointMap` and (un)differentiated
        // `target`, we can also introduce a new `generator` that is
        // differentiated with respect to `s`.
        let diff_generator = self.generator.differentiate(s).map_err(|e| {
            generic_expression_error(
                "ExpAdjointMap::differentiate() failed for generator",
                self,
                Some(Box::new(e)),
            )
        })?;

        let new_deriv = self.derivative.with_added_perturbation(s);

        if is_expr_type::<ZeroOperator>(&diff_generator) {
            return if is_expr_type::<ZeroOperator>(&diff_result) {
                Ok(diff_result)
            } else {
                self.with_result_and_derivative(diff_result, new_deriv).build()
            };
        }

        let mut terms = Vec::new();
        if !is_expr_type::<ZeroOperator>(&diff_result) {
            terms.push(diff_result);
        }

        // `adj_maps` contains `AdjointMap`'s that should be extracted from the
        // exponential adjoint map due to maximum commutator order
        let mut adj_maps = Vec::new();

        let adjoint_mode = if self.generator_derivative_commute {
            Some(AdjointMode::Commutative)
        } else {
            Some(AdjointMode::Symmetric)
        };

        if let Some(mat_add) = downcast_from_arc::<MatrixAdd>(&self.result) {
            for term in mat_add.terms() {
                if let Some(adj_map) = downcast_from_arc::<AdjointMap>(term) {
                    // Check commutator order
                    if adj_map.generators().len() as u32 + 1 < self.max_commutator_order {
                        terms.push(
                            adj_map.with_added_generator(diff_generator.clone(), adjoint_mode)?,
                        );
                    } else {
                        adj_maps.push(term.clone());
                    }
                } else {
                    // `term` is a differentiated `target`
                    terms.push(AdjointMap::new(
                        vec![diff_generator.clone()],
                        term.clone(),
                        Some(self.left_action),
                        adjoint_mode,
                    )?);
                }
            }
        } else if let Some(adj_map) = downcast_from_arc::<AdjointMap>(&self.result) {
            if adj_map.generators().len() as u32 + 1 < self.max_commutator_order {
                terms.push(adj_map.with_added_generator(diff_generator.clone(), adjoint_mode)?);
            } else {
                adj_maps.push(self.result.clone());
            }
        } else {
            // `result` is the undifferentiated `target`
            terms.push(AdjointMap::new(
                vec![diff_generator],
                self.result.clone(),
                Some(self.left_action),
                adjoint_mode,
            )?);
        }

        let new_ead_map =
            self.with_result_and_derivative(MatrixAdd::new(terms)?, new_deriv).build()?;

        if adj_maps.is_empty() {
            Ok(new_ead_map)
        } else {
            adj_maps.push(new_ead_map);
            MatrixAdd::new(adj_maps)
        }
    }
}

impl PartialEq for ExpAdjointMap {
    fn eq(&self, other: &Self) -> bool {
        // We also compare `result`, which may change after `substitute_zero_perturbations()`
        &self.generator == &other.generator
            && self.generator_derivative_commute == other.generator_derivative_commute
            && &self.target == &other.target
            && self.is_time_evolution == other.is_time_evolution
            && self.left_action == other.left_action
            && self.max_commutator_order == other.max_commutator_order
            && self.derivative == other.derivative
            && &self.result == &other.result
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
