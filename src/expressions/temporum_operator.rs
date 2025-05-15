use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{MatrixMul, OneElecOperator, WfnParameter, ZeroOperator};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{
    NumberTolerance, downcast_from_arc, downcast_from_ref, expression_error,
    generic_expression_error, is_expr_type, is_zero_expr, negate_expr, sum_pert_frequencies,
    unreachable_error,
};

/// A TemporumOperator represents i*d/dt (forward) or -i*d/dt (backward) acting
/// on an `argument`.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct TemporumOperator {
    is_forward: bool,
    argument: Arc<dyn Expr>,
}

impl TemporumOperator {
    #[inline]
    pub fn builder(argument: Arc<dyn Expr>) -> TemporumOperatorBuilder {
        TemporumOperatorBuilder {
            is_forward: true,
            argument,
        }
    }

    #[inline]
    fn builder_from(&self, argument: Arc<dyn Expr>) -> TemporumOperatorBuilder {
        TemporumOperatorBuilder {
            is_forward: self.is_forward,
            argument,
        }
    }

    #[inline]
    pub fn is_forward(&self) -> bool {
        self.is_forward
    }

    #[inline]
    pub fn argument(&self) -> &Arc<dyn Expr> {
        &self.argument
    }

    #[inline]
    pub fn derivative(&self) -> Result<&PertMultichain, TinnedError> {
        if let Some(op) = downcast_from_arc::<OneElecOperator>(&self.argument) {
            Ok(op.derivative())
        } else if let Some(wfn) = downcast_from_arc::<WfnParameter>(&self.argument) {
            Ok(wfn.derivative())
        } else {
            Err(unreachable_error(
                "TemporumOperator::frequency() gets an argument neither OneElecOperator nor WfnParameter",
                &self.argument,
                None,
            ))
        }
    }

    // For unperturbed `argument`, the function `frequency()` should return
    // zero number
    #[inline]
    pub fn frequency(&self) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.is_forward {
            sum_pert_frequencies(self.derivative()?)
        } else {
            negate_expr(sum_pert_frequencies(self.derivative()?)?)
        }
    }
}

#[derive(Debug)]
pub struct TemporumOperatorBuilder {
    is_forward: bool,
    argument: Arc<dyn Expr>,
}

impl TemporumOperatorBuilder {
    #[inline]
    pub fn is_forward(mut self, is_forward: bool) -> Self {
        self.is_forward = is_forward;
        self
    }

    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.argument.is_scalar() {
            return Err(expression_error(
                "TemporumOperatorBuilder::build() - scalar argument",
                &self.argument,
                None,
            ));
        }

        if is_expr_type::<ZeroOperator>(&self.argument) {
            return Ok(self.argument);
        }

        if is_expr_type::<OneElecOperator>(&self.argument)
            || is_expr_type::<WfnParameter>(&self.argument)
        {
            Ok(crate::internal::intern_expr(Arc::new(TemporumOperator {
                is_forward: self.is_forward,
                argument: self.argument,
            })))
        } else {
            Err(expression_error(
                "TemporumOperatorBuilder::build() - unsupported argument type",
                &self.argument,
                None,
            ))
        }
    }
}

impl ExprInternal for TemporumOperator {
    impl_expr_internal_methods!(TemporumOperator);

    #[inline]
    fn hash_key(&self) -> String {
        format!("TemporumOperator({}; {})", self.is_forward, self.argument.hash_key())
    }

    #[inline]
    fn total_order(&self) -> u32 {
        self.argument.total_order()
    }

    #[inline]
    fn match_for_find_all(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<TemporumOperator>(other) {
            self.argument.match_for_find_all(&op.argument)
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for TemporumOperator {
    impl_unary_expr_common_methods!(
        TemporumOperator,
        argument,
        False,
        false,
        |this: &TemporumOperator, arg| this.builder_from(arg).build()
    );

    #[inline]
    fn clean_temporum(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let frequency = self.frequency()?;

        if is_zero_expr(&frequency, freq_tol) {
            Ok(ZeroOperator::new())
        } else {
            MatrixMul::new(vec![frequency, self.argument.clone()])
        }
    }

    #[inline]
    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_arg = self.argument.differentiate(s).map_err(|e| {
            generic_expression_error(
                "TemporumOperator::differentiate() failed for argument",
                self,
                Some(Box::new(e)),
            )
        })?;

        self.builder_from(diff_arg).build()
    }
}

impl PartialEq for TemporumOperator {
    fn eq(&self, other: &Self) -> bool {
        self.is_forward == other.is_forward && &self.argument == &other.argument
    }
}

impl Eq for TemporumOperator {}

impl std::fmt::Display for TemporumOperator {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}({})",
            if self.is_forward {
                "i*d/dt"
            } else {
                "-i*d/dt"
            },
            self.argument,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::one_elec_operator::test_utils::make_one_elec_operator;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::is_one_expr;

    test_struct_safety!(TemporumOperator);

    test_thread_interning!(
        TemporumOperator::builder(make_one_elec_operator("1el")).build().unwrap()
    );

    #[test]
    fn test_impl_expr() {
        let is_forward = true;
        let argument = make_one_elec_operator("");
        let op1 =
            TemporumOperator::builder(argument.clone()).is_forward(is_forward).build().unwrap();

        let op = downcast_from_arc::<TemporumOperator>(&op1).unwrap();
        assert_eq!(
            op,
            &TemporumOperator {
                is_forward,
                argument: argument.clone()
            }
        );

        assert_eq!(op.is_forward(), is_forward);
        assert_eq!(op.argument(), &argument.clone());

        let op2 = op.builder_from(argument.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        assert_eq!(
            op1.hash_key(),
            format!("TemporumOperator({}; {})", is_forward, argument.hash_key())
        );
        assert!(!op1.is_scalar());
        assert_eq!(
            format!("{}", op1),
            format!(
                "{}({})",
                if is_forward {
                    "i*d/dt"
                } else {
                    "-i*d/dt"
                },
                argument,
            )
        );

        let op3 =
            TemporumOperator::builder(argument.clone()).is_forward(!is_forward).build().unwrap();
        let op4 = TemporumOperator::builder(make_wfn_parameter("")).build().unwrap();

        assert_ne!(&op1, &op3);
        assert_ne!(&op1, &op4);
    }

    #[test]
    fn test_differentiation() {
        let is_forward = true;
        let mut argument = make_one_elec_operator("");
        let op1 =
            TemporumOperator::builder(argument.clone()).is_forward(is_forward).build().unwrap();

        let p = make_perturbation_symbol(4u32, 4u32);
        let diff_op1 = op1.differentiate(&p).unwrap();
        let diff_arg = argument.differentiate(&p).unwrap();

        if is_expr_type::<ZeroOperator>(&diff_arg) {
            assert!(is_expr_type::<ZeroOperator>(&diff_op1));
        } else {
            let diff_cast = downcast_from_arc::<TemporumOperator>(&diff_op1).unwrap();

            assert_eq!(diff_cast.argument(), &diff_arg);
        }

        argument = make_wfn_parameter("");
        let op2 = TemporumOperator::builder(argument.clone()).build().unwrap();
        let diff_op2 = op2.differentiate(&p).unwrap();
        let diff_cast = downcast_from_arc::<TemporumOperator>(&diff_op2).unwrap();

        assert_eq!(diff_cast.argument(), &argument.differentiate(&p).unwrap());
    }

    #[test]
    fn test_serialization() {
        let mut op = TemporumOperator::builder(make_one_elec_operator("")).build().unwrap();
        let mut json = serde_json::to_string(&op).unwrap();
        let mut deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);

        op = TemporumOperator::builder(make_wfn_parameter("")).build().unwrap();
        json = serde_json::to_string(&op).unwrap();
        deserialized = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let op1 = TemporumOperator::builder(make_one_elec_operator("1el")).build().unwrap();

        assert!(is_expr_type::<TemporumOperator>(&op1));
        assert!(!is_zero_expr(&op1, None));
        assert!(!is_one_expr(&op1, None));

        let op2 = TemporumOperator::builder(make_one_elec_operator("1el")).build().unwrap();
        let op3 = TemporumOperator::builder(make_one_elec_operator("")).build().unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));
    }
}
