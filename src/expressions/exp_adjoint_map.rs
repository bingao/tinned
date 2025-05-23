use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{AdjointMap, MatrixAdd, ZeroOperator};
use crate::internal::intern_expr;
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{
    NumberTolerance, differentiate_expr, downcast_from_arc, downcast_from_ref, expression_error,
    generic_expression_error, is_expr_type, unreachable_error,
};

// Exponential adjoint map (or conjugation operation in Lie algebra):
// exp(ad_{X})(Y) = exp(X)*Y*exp(-X) (`left_action` is `true`), or
// exp(ad_{-X})(Y) = exp(-X)*Y*exp(X) (`left_action` is `false`).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct ExpAdjointMap {
    generator: Arc<dyn Expr>,
    target: Arc<dyn Expr>,
    left_action: bool,
    max_fold: u32,
    is_zero_strength: bool,
    // `result` contains differentiated expression of exponential adjoint map
    result: Arc<dyn Expr>,
    derivative: PertMultichain,
}

impl ExpAdjointMap {
    #[inline]
    pub fn builder(generator: Arc<dyn Expr>, target: Arc<dyn Expr>) -> ExpAdjointMapBuilder {
        ExpAdjointMapBuilder {
            generator,
            target,
            left_action: None,
            max_fold: None,
            is_zero_strength: Some(false),
            result: None,
            derivative: None,
        }
    }

    #[inline]
    fn with_result(
        &self,
        result: Arc<dyn Expr>,
        is_zero_strength: Option<bool>,
    ) -> ExpAdjointMapBuilder {
        ExpAdjointMapBuilder {
            generator: self.generator.clone(),
            target: self.target.clone(),
            left_action: Some(self.left_action),
            max_fold: Some(self.max_fold),
            is_zero_strength,
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
            target: self.target.clone(),
            left_action: Some(self.left_action),
            max_fold: Some(self.max_fold),
            is_zero_strength: Some(self.is_zero_strength),
            result: Some(result),
            derivative: Some(derivative),
        }
    }

    #[inline]
    pub fn generator(&self) -> &Arc<dyn Expr> {
        &self.generator
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
    pub fn max_fold(&self) -> u32 {
        self.max_fold
    }

    #[inline]
    pub fn is_zero_strength(&self) -> bool {
        self.is_zero_strength
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
    target: Arc<dyn Expr>,
    left_action: Option<bool>,
    max_fold: Option<u32>,
    is_zero_strength: Option<bool>,
    result: Option<Arc<dyn Expr>>,
    derivative: Option<PertMultichain>,
}

impl ExpAdjointMapBuilder {
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
    //fn is_zero_strength(mut self, is_zero_strength: bool) -> Self {
    //    self.is_zero_strength = Some(is_zero_strength);
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

        let left_action = self.left_action.unwrap_or(true);
        let max_fold = self.max_fold.unwrap_or(u32::MAX);
        let is_zero_strength = self.is_zero_strength.unwrap_or(false);

        // Undifferentiated expression of exponential adjoint map is simply `target`
        let result = self.result.unwrap_or(self.target.clone());
        let derivative = self.derivative.unwrap_or(PertMultichain::new());

        Ok(intern_expr(Arc::new(ExpAdjointMap {
            generator: self.generator,
            target: self.target,
            left_action,
            max_fold,
            is_zero_strength,
            result,
            derivative,
        })))
    }
}

impl ExprInternal for ExpAdjointMap {
    impl_unary_expr_internal_methods!(ExpAdjointMap, result, true, |this: &ExpAdjointMap, arg| {
        this.with_result(arg, Some(this.is_zero_strength)).build()
    });

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "ExpAdjointMap({}; {}; {}; {}; {}; {}; [{}])",
            self.left_action,
            self.max_fold,
            self.is_zero_strength,
            self.generator.hash_key(),
            self.target.hash_key(),
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
            // We find exponential adjoint maps with `left_action` and
            // `is_zero_strength` either `true` or `false`
            self.max_fold == op.max_fold
                && self.generator.deep_eq_superchains(&op.generator)
                && self.target.deep_eq_superchains(&op.target)
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
                && self.is_zero_strength == op.is_zero_strength
                && &self.generator == &op.generator
                && &self.target == &op.target
                && self.derivative.is_subchain(&op.derivative)
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for ExpAdjointMap {
    impl_unary_expr_common_methods!(ExpAdjointMap, result, False, |this: &ExpAdjointMap, arg| this
        .with_result(arg, Some(this.is_zero_strength))
        .build());

    #[inline]
    fn clean_temporum(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.is_zero_strength {
            return Ok(self.clone_expr());
        }

        let result = self.result.clean_temporum(freq_tol.clone()).map_err(|e| {
            generic_expression_error(
                "ExpAdjointMap::clean_temporum() failed for result",
                self,
                Some(Box::new(e)),
            )
        })?;

        self.with_result(result, Some(true)).build()
    }

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

        let mut ad_maps = Vec::new();

        if let Some(mat_add) = downcast_from_arc::<MatrixAdd>(&self.result) {
            for term in mat_add.terms() {
                if let Some(ad_map) = downcast_from_arc::<AdjointMap>(term) {
                    // Check folds of commutators
                    if ad_map.generators().len() as u32 + 1 < self.max_fold {
                        terms.push(ad_map.with_added_generator(diff_generator.clone())?);
                    } else {
                        ad_maps.push(term.clone());
                    }
                } else {
                    // differentiated `target`
                    terms.push(AdjointMap::new(
                        vec![diff_generator.clone()],
                        term.clone(),
                        Some(self.left_action),
                    )?);
                }
            }
        } else {
            // `result` is undifferentiated `target`
            terms.push(AdjointMap::new(
                vec![diff_generator],
                self.result.clone(),
                Some(self.left_action),
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
        // We also compare `result`, which may change after `clean_temporum()`
        &self.generator == &other.generator
            && &self.target == &other.target
            && self.left_action == other.left_action
            && self.max_fold == other.max_fold
            && self.is_zero_strength == other.is_zero_strength
            && self.derivative == other.derivative
            && &self.result == &other.result
    }
}

impl Eq for ExpAdjointMap {}

impl std::fmt::Display for ExpAdjointMap {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.left_action {
            write!(f, "exp(ad[{}", self.generator)?;
        } else {
            write!(f, "exp(ad[-({})", self.generator)?;
        }

        if self.max_fold < u32::MAX {
            write!(f, "; {}", self.max_fold)?;
        }

        write!(f, "])({}; {})^{}", self.target, self.is_zero_strength, self.derivative)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    test_struct_safety!(ExpAdjointMap);
}
