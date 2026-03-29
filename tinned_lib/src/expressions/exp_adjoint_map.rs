use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{AdjointMap, AdjointMode, MatrixAdd, TemporumOperator, ZeroOperator};
use crate::internal::intern_expr;
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{
    NumberTolerance, downcast_from_arc, expression_error, generic_expression_error, is_expr_type,
};

// Exponential adjoint map (or conjugation operation in Lie algebra):
// exp(ad_{X})(Y) = exp(X)*Y*exp(-X) (`left_action` is `true`), or
// exp(ad_{-X})(Y) = exp(-X)*Y*exp(X) (`left_action` is `false`).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct ExpAdjointMap {
    generator: Arc<dyn Expr>,
    // Whether generator and its derivatives commute
    generator_derivative_commute: bool,
    target: Arc<dyn Expr>,
    is_temporum: bool,
    left_action: bool,
    max_fold: u32,
    zero_rules_applied: bool,
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
            is_temporum: false,
            left_action: None,
            max_fold: None,
            zero_rules_applied: Some(false),
            result: None,
            derivative: None,
        }
    }

    #[inline]
    pub fn builder_temporum(
        generator: Arc<dyn Expr>,
        is_forward: bool,
        generator_derivative_commute: Option<bool>,
    ) -> ExpAdjointMapBuilder {
        let target =
            match TemporumOperator::builder(generator.clone()).is_forward(is_forward).build() {
                Ok(e) => e,
                Err(e) => panic!("ExpAdjointMap::builder_temporum() encounters: {e}"),
            };

        ExpAdjointMapBuilder {
            generator,
            generator_derivative_commute,
            target,
            is_temporum: true,
            left_action: None,
            max_fold: None,
            zero_rules_applied: Some(false),
            result: None,
            derivative: None,
        }
    }

    // This function and `with_result_and_derivative()` are only used inside this file
    #[inline]
    fn with_result(
        &self,
        result: Arc<dyn Expr>,
        zero_rules_applied: Option<bool>,
    ) -> ExpAdjointMapBuilder {
        ExpAdjointMapBuilder {
            generator: self.generator.clone(),
            generator_derivative_commute: Some(self.generator_derivative_commute),
            target: self.target.clone(),
            is_temporum: self.is_temporum,
            left_action: Some(self.left_action),
            max_fold: Some(self.max_fold),
            zero_rules_applied,
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
            is_temporum: self.is_temporum,
            left_action: Some(self.left_action),
            max_fold: Some(self.max_fold),
            zero_rules_applied: Some(self.zero_rules_applied),
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
    pub fn is_temporum(&self) -> bool {
        self.is_temporum
    }

    #[inline]
    pub fn left_action(&self) -> bool {
        self.left_action
    }

    #[inline]
    pub fn max_fold(&self) -> u32 {
        self.max_fold
    }

    #[inline]
    pub fn zero_rules_applied(&self) -> bool {
        self.zero_rules_applied
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
    is_temporum: bool,
    left_action: Option<bool>,
    max_fold: Option<u32>,
    zero_rules_applied: Option<bool>,
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
    pub fn max_fold(mut self, max_fold: u32) -> Self {
        self.max_fold = Some(max_fold);
        self
    }

    //#[inline]
    //fn zero_rules_applied(mut self, zero_rules_applied: bool) -> Self {
    //    self.zero_rules_applied = Some(zero_rules_applied);
    //    self
    //}

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
        let max_fold = self.max_fold.unwrap_or(u32::MAX);
        let zero_rules_applied = self.zero_rules_applied.unwrap_or(false);

        // Undifferentiated expression of exponential adjoint map is simply `target`
        let result = self.result.unwrap_or(self.target.clone());
        let derivative = self.derivative.unwrap_or(PertMultichain::new());

        Ok(intern_expr(Arc::new(ExpAdjointMap {
            generator: self.generator,
            generator_derivative_commute,
            target: self.target,
            is_temporum: self.is_temporum,
            left_action,
            max_fold,
            zero_rules_applied,
            result,
            derivative,
        })))
    }
}

impl ExprInternal for ExpAdjointMap {
    // It is more appropriate to set `zero_rules_applied` as false after the
    // functions `replace_expr_children`, `retain_expr_fields` and `replace_expr_self`
    impl_unary_expr_internal_methods!(
        ExpAdjointMap,
        False,
        result,
        true,
        |this: &ExpAdjointMap, arg| { this.with_result(arg, Some(false)).build() }
    );

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "ExpAdjointMap({}; {}; {}; {}; {}; {}; {}; {}; [{}])",
            self.left_action,
            self.max_fold,
            self.zero_rules_applied,
            self.generator.hash_key(),
            self.generator_derivative_commute,
            self.target.hash_key(),
            self.is_temporum,
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
            // We treat exponential adjoint maps with different `left_action`
            // and `zero_rules_applied` equally
            self.max_fold == op.max_fold
                && self.generator.deep_eq_superchains(&op.generator)
                && self.generator_derivative_commute == op.generator_derivative_commute
                && self.target.deep_eq_superchains(&op.target)
                && self.is_temporum == op.is_temporum
                && self.derivative.is_subchain(&op.derivative)
        } else {
            false
        }
    }

    #[inline]
    fn eq_by_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<ExpAdjointMap>(other) {
            self.left_action == op.left_action
                && self.max_fold == op.max_fold
                && self.zero_rules_applied == op.zero_rules_applied
                && &self.generator == &op.generator
                && self.generator_derivative_commute == op.generator_derivative_commute
                && &self.target == &op.target
                && self.is_temporum == op.is_temporum
                && self.derivative.is_subchain(&op.derivative)
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for ExpAdjointMap {
    impl_unary_expr_common_methods!(ExpAdjointMap, False, result, |this: &ExpAdjointMap, arg| this
        .with_result(arg, Some(this.zero_rules_applied))
        .build());

    #[inline]
    fn apply_zero_rules(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.zero_rules_applied {
            return Ok(self.clone_expr());
        }

        let result = self.result.apply_zero_rules(freq_tol.clone()).map_err(|e| {
            generic_expression_error(
                "ExpAdjointMap::apply_zero_rules() failed for result",
                self,
                Some(Box::new(e)),
            )
        })?;

        self.with_result(result, Some(true)).build()
    }

    // Equation (47)
    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        // `result` is (i) an `MatrixAdd` of `AdjointMap`'s and differentiated
        // `target`, or (ii) undifferentiated `target`.  We first differentiate
        // `result` with respect to `s`, which gives us all differentiated
        // terms (including the differentiation on each `AdjointMap`'s field
        // `target`) with number of generators fixed for each `AdjointMap`.
        let diff_result = self.result.differentiate(s).map_err(|e| {
            generic_expression_error(
                "ExpAdjointMap::differentiate() failed for result",
                self,
                Some(Box::new(e)),
            )
        })?;

        // The first order derivative must be performed on the temporum
        // operator +/-i*d/dt
        if self.is_temporum && self.derivative.is_empty() {
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

        // `ad_maps` contains `AdjointMap`'s that should be extracted from the
        // exponential adjoint map due to maximum folds of commutators
        let mut ad_maps = Vec::new();

        let adjoint_mode = if self.generator_derivative_commute {
            Some(AdjointMode::Commutative)
        } else {
            Some(AdjointMode::Symmetric)
        };

        if let Some(mat_add) = downcast_from_arc::<MatrixAdd>(&self.result) {
            for term in mat_add.terms() {
                if let Some(ad_map) = downcast_from_arc::<AdjointMap>(term) {
                    // Check folds of commutators
                    if ad_map.generators().len() as u32 + 1 < self.max_fold {
                        terms.push(
                            ad_map.with_added_generator(diff_generator.clone(), adjoint_mode)?,
                        );
                    } else {
                        ad_maps.push(term.clone());
                    }
                } else {
                    // differentiated `target`
                    terms.push(AdjointMap::new(
                        vec![diff_generator.clone()],
                        term.clone(),
                        Some(self.left_action),
                        adjoint_mode,
                    )?);
                }
            }
        } else {
            // `result` is undifferentiated `target`
            terms.push(AdjointMap::new(
                vec![diff_generator],
                self.result.clone(),
                Some(self.left_action),
                adjoint_mode,
            )?);
        }

        let new_ead_map =
            self.with_result_and_derivative(MatrixAdd::new(terms)?, new_deriv).build()?;

        if ad_maps.is_empty() {
            Ok(new_ead_map)
        } else {
            ad_maps.push(new_ead_map);
            MatrixAdd::new(ad_maps)
        }
    }
}

impl PartialEq for ExpAdjointMap {
    fn eq(&self, other: &Self) -> bool {
        // We also compare `result`, which may change after `apply_zero_rules()`
        &self.generator == &other.generator
            && self.generator_derivative_commute == other.generator_derivative_commute
            && &self.target == &other.target
            && self.is_temporum == other.is_temporum
            && self.left_action == other.left_action
            && self.max_fold == other.max_fold
            && self.zero_rules_applied == other.zero_rules_applied
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

        if self.max_fold < u32::MAX {
            write!(f, "; {}", self.max_fold)?;
        }

        write!(f, "])({}; {})^{}", self.target, self.zero_rules_applied, self.derivative)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    test_struct_safety!(ExpAdjointMap);
}
