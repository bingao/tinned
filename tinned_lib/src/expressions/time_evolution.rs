use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{MatrixMul, OneElecMatrix, ResidueParameter, WfnParameter, ZeroOperator};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{
    NumberTolerance, downcast_from_arc, expression_error, generic_expression_error, is_expr_type,
    is_zero_expr, negate_expr, sum_pert_frequencies, unreachable_error,
};

/// A TimeEvolution represents i*d/dt (forward) or -i*d/dt (backward) acting
/// on an `argument`.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct TimeEvolution {
    is_forward: bool,
    argument: Arc<dyn Expr>,
}

impl TimeEvolution {
    #[inline]
    pub fn builder(argument: Arc<dyn Expr>) -> TimeEvolutionBuilder {
        TimeEvolutionBuilder {
            is_forward: true,
            argument,
        }
    }

    #[inline]
    fn with_argument(&self, argument: Arc<dyn Expr>) -> TimeEvolutionBuilder {
        TimeEvolutionBuilder {
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
        if let Some(op) = downcast_from_arc::<OneElecMatrix>(&self.argument) {
            Ok(op.derivative())
        } else if let Some(wfn) = downcast_from_arc::<WfnParameter>(&self.argument) {
            Ok(wfn.derivative())
        } else if let Some(residue) = downcast_from_arc::<ResidueParameter>(&self.argument) {
            residue.derivative()
        } else {
            Err(unreachable_error(
                "TimeEvolution::derivative() gets an argument neither OneElecMatrix nor WfnParameter",
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
pub struct TimeEvolutionBuilder {
    is_forward: bool,
    argument: Arc<dyn Expr>,
}

impl TimeEvolutionBuilder {
    #[inline]
    pub fn is_forward(mut self, is_forward: bool) -> Self {
        self.is_forward = is_forward;
        self
    }

    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.argument.is_scalar() {
            return Err(expression_error(
                "TimeEvolutionBuilder::build() - scalar argument",
                &self.argument,
                None,
            ));
        }

        if is_expr_type::<ZeroOperator>(&self.argument) {
            return Ok(self.argument);
        }

        if is_expr_type::<OneElecMatrix>(&self.argument)
            || is_expr_type::<WfnParameter>(&self.argument)
            || is_expr_type::<ResidueParameter>(&self.argument)
        {
            Ok(crate::internal::intern_expr(Arc::new(TimeEvolution {
                is_forward: self.is_forward,
                argument: self.argument,
            })))
        } else {
            Err(expression_error(
                "TimeEvolutionBuilder::build() - unsupported argument type",
                &self.argument,
                None,
            ))
        }
    }
}

impl ExprInternal for TimeEvolution {
    impl_unary_expr_internal_methods!(
        TimeEvolution,
        False,
        argument,
        false,
        |this: &TimeEvolution, arg| this.with_argument(arg).build()
    );

    #[inline]
    fn hash_key(&self) -> String {
        format!("TimeEvolution({}; {})", self.is_forward, self.argument.hash_key())
    }

    #[inline]
    fn total_order(&self) -> u32 {
        self.argument.total_order()
    }

    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<TimeEvolution>(other) {
            self.argument.deep_eq_superchains(&op.argument)
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for TimeEvolution {
    impl_unary_expr_common_methods!(TimeEvolution, False, argument, |this: &TimeEvolution, arg| {
        this.with_argument(arg).build()
    });

    #[inline]
    fn has_unperturbed_term(&self) -> bool {
        false
    }

    #[inline]
    fn substitute_zero_perturbations(
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
    fn differentiate(&self, s: Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_arg = self.argument.differentiate(s).map_err(|e| {
            generic_expression_error(
                "TimeEvolution::differentiate() failed for argument",
                self,
                Some(Box::new(e)),
            )
        })?;

        self.with_argument(diff_arg).build()
    }
}

impl PartialEq for TimeEvolution {
    fn eq(&self, other: &Self) -> bool {
        self.is_forward == other.is_forward && &self.argument == &other.argument
    }
}

impl Eq for TimeEvolution {}

impl std::fmt::Display for TimeEvolution {
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
    use crate::expressions::one_elec_matrix::test_utils::make_one_elec_matrix;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::is_one_expr;

    test_struct_safety!(TimeEvolution);

    test_thread_interning!(
        TimeEvolution::builder(make_one_elec_matrix("1el", false)).build().unwrap()
    );

    #[test]
    fn test_impl_expr() {
        let is_forward = true;
        let argument = make_one_elec_matrix("", false);
        let op1 = TimeEvolution::builder(argument.clone()).is_forward(is_forward).build().unwrap();

        let op = downcast_from_arc::<TimeEvolution>(&op1).unwrap();
        assert_eq!(
            op,
            &TimeEvolution {
                is_forward,
                argument: argument.clone()
            }
        );

        assert_eq!(op.is_forward(), is_forward);
        assert_eq!(op.argument(), &argument.clone());

        let op2 = op.with_argument(argument.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        assert_eq!(
            op1.hash_key(),
            format!("TimeEvolution({}; {})", is_forward, argument.hash_key())
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

        let op3 = TimeEvolution::builder(argument.clone()).is_forward(!is_forward).build().unwrap();
        let op4 = TimeEvolution::builder(make_wfn_parameter("")).build().unwrap();

        assert_ne!(&op1, &op3);
        assert_ne!(&op1, &op4);
    }

    #[test]
    fn test_differentiation() {
        let is_forward = true;
        let mut argument = make_one_elec_matrix("", false);
        let op1 = TimeEvolution::builder(argument.clone()).is_forward(is_forward).build().unwrap();

        let p = make_perturbation_symbol(4u32, 4u32);
        let diff_op1 = op1.differentiate(p.clone()).unwrap();
        let diff_arg = argument.differentiate(p.clone()).unwrap();

        if is_expr_type::<ZeroOperator>(&diff_arg) {
            assert!(is_expr_type::<ZeroOperator>(&diff_op1));
        } else {
            let diff_cast = downcast_from_arc::<TimeEvolution>(&diff_op1).unwrap();

            assert_eq!(diff_cast.argument(), &diff_arg);
        }

        argument = make_wfn_parameter("");
        let op2 = TimeEvolution::builder(argument.clone()).build().unwrap();
        let diff_op2 = op2.differentiate(p.clone()).unwrap();
        let diff_cast = downcast_from_arc::<TimeEvolution>(&diff_op2).unwrap();

        assert_eq!(diff_cast.argument(), &argument.differentiate(p).unwrap());
    }

    #[test]
    fn test_serialization() {
        let mut op = TimeEvolution::builder(make_one_elec_matrix("", false)).build().unwrap();
        let mut json = serde_json::to_string(&op).unwrap();
        let mut deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);

        op = TimeEvolution::builder(make_wfn_parameter("")).build().unwrap();
        json = serde_json::to_string(&op).unwrap();
        deserialized = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let op1 = TimeEvolution::builder(make_one_elec_matrix("1el", false)).build().unwrap();

        assert!(is_expr_type::<TimeEvolution>(&op1));
        assert!(!is_zero_expr(&op1, None));
        assert!(!is_one_expr(&op1, None));

        let op2 = TimeEvolution::builder(make_one_elec_matrix("1el", false)).build().unwrap();
        let op3 = TimeEvolution::builder(make_one_elec_matrix("", false)).build().unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));
    }
}
