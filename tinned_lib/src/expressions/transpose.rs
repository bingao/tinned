use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{Conjugate, MatrixMul, ZeroOperator};
use crate::internal::{intern_expr, transform_unary_any_zero};
use crate::perturbations::Perturbation;
use crate::public::{
    NumberTolerance, downcast_from_arc, expression_error, is_expr_type, is_one_expr,
};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Transpose {
    argument: Arc<dyn Expr>,
    is_hermitian: bool,
}

impl Transpose {
    pub fn new(argument: Arc<dyn Expr>, is_hermitian: bool) -> Result<Arc<dyn Expr>, TinnedError> {
        if argument.is_scalar() {
            return Err(expression_error(
                "Transpose::new() got a scalar argument",
                &argument,
                None,
            ));
        }

        if is_expr_type::<ZeroOperator>(&argument) {
            return Ok(argument);
        } else if let Some(conj) = downcast_from_arc::<Conjugate>(&argument) {
            return Transpose::new(conj.argument().clone(), !is_hermitian);
        } else if let Some(trans) = downcast_from_arc::<Transpose>(&argument) {
            return if is_hermitian == trans.is_hermitian {
                Ok(trans.argument().clone())
            } else {
                Conjugate::new(trans.argument().clone())
            };
        } else if let Some(mat_mul) = downcast_from_arc::<MatrixMul>(&argument) {
            if is_one_expr(mat_mul.coefficient(), None) {
                return Ok(intern_expr(Arc::new(Self {
                    argument,
                    is_hermitian,
                })));
            }

            let new_coefficient = if is_hermitian {
                Conjugate::new(mat_mul.coefficient().clone())?
            } else {
                mat_mul.coefficient().clone()
            };
            let new_argument = MatrixMul::new(mat_mul.factors().to_vec())?;

            return MatrixMul::new(vec![
                new_coefficient,
                intern_expr(Arc::new(Self {
                    argument: new_argument,
                    is_hermitian,
                })),
            ]);
        }

        Ok(intern_expr(Arc::new(Self {
            argument,
            is_hermitian,
        })))
    }

    #[inline]
    pub fn argument(&self) -> &Arc<dyn Expr> {
        &self.argument
    }

    #[inline]
    pub fn is_hermitian(&self) -> bool {
        self.is_hermitian
    }
}

impl ExprInternal for Transpose {
    impl_unary_expr_internal_methods!(
        Transpose,
        false,
        argument,
        false,
        |this: &Transpose, arg| Self::new(arg, this.is_hermitian)
    );

    #[inline]
    fn hash_key(&self) -> String {
        format!("Transpose({}; {})", self.argument.hash_key(), self.is_hermitian)
    }

    #[inline]
    fn expr_order(&self) -> u32 {
        self.argument.expr_order()
    }

    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(expr) = downcast_from_arc::<Transpose>(other) {
            self.argument.deep_eq_superchains(&expr.argument)
        } else {
            false
        }
    }

    // Like most unary `Expr` types, `Transpose` itself does not hold any
    // derivative. So we use equality comparison on the whole expression, i.e.
    // we do not override the method `eq_by_superchains()` of the trait
    // `ExprInternal`.
}

#[typetag::serde]
impl Expr for Transpose {
    impl_unary_expr_common_methods!(Transpose, false, argument, |this: &Transpose, arg| {
        Self::new(arg, this.is_hermitian)
    });

    #[inline]
    fn has_unperturbed_term(&self) -> bool {
        self.argument.has_unperturbed_term()
    }

    #[inline]
    fn substitute_zero_perturbations(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        transform_unary_any_zero(
            self,
            &self.argument,
            |arg: &Arc<dyn Expr>| arg.substitute_zero_perturbations(freq_tol),
            "Transpose::substitute_zero_perturbations() failed for argument",
            |arg| Self::new(arg, self.is_hermitian),
            || Ok(ZeroOperator::new()),
        )
    }

    #[inline]
    fn differentiate(&self, s: Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        transform_unary_any_zero(
            self,
            &self.argument,
            |arg: &Arc<dyn Expr>| arg.differentiate(s),
            "Transpose::differentiate() failed for argument",
            |arg| Self::new(arg, self.is_hermitian),
            || Ok(ZeroOperator::new()),
        )
    }
}

impl PartialEq for Transpose {
    fn eq(&self, other: &Self) -> bool {
        &self.argument == &other.argument && self.is_hermitian == other.is_hermitian
    }
}

impl Eq for Transpose {}

impl std::fmt::Display for Transpose {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_hermitian {
            write!(f, "{}^H", self.argument)
        } else {
            write!(f, "{}^T", self.argument)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::ao_two_elec_matrix::test_utils::make_ao_two_elec_matrix;
    use crate::expressions::symbol::test_utils::make_symbol;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::public::is_zero_expr;

    test_unary_oper_properties!(Transpose, |arg| Transpose::new(arg, false));

    #[test]
    fn test_impl_expr() {
        let standard_transpose = false;
        let hermitian_transpose = true;

        let op0 = Transpose::new(ZeroOperator::new(), standard_transpose).unwrap();
        assert!(is_zero_expr(&op0, None));

        let arg_2el = make_ao_two_elec_matrix("", None);
        let op1 = Transpose::new(arg_2el.clone(), standard_transpose).unwrap();

        let op = downcast_from_arc::<Transpose>(&op1).unwrap();
        assert_eq!(
            op,
            &Transpose {
                argument: arg_2el.clone(),
                is_hermitian: standard_transpose,
            }
        );
        assert_eq!(op.argument(), &arg_2el);

        let op2 = Transpose::new(arg_2el.clone(), standard_transpose).unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        assert_eq!(
            op1.hash_key(),
            format!("Transpose({}; {})", arg_2el.hash_key(), standard_transpose)
        );
        assert!(!op1.is_scalar());
        assert_eq!(format!("{}", op1), format!("{}^T", arg_2el));

        let mut argument = Conjugate::new(arg_2el.clone()).unwrap();
        let op3 = Transpose::new(argument, standard_transpose).unwrap();
        assert_ne!(&op1, &op3);
        assert_eq!(&op3, &Transpose::new(arg_2el.clone(), hermitian_transpose).unwrap());

        let op4 = Transpose::new(op2, standard_transpose).unwrap();
        assert_ne!(&op1, &op4);
        assert_eq!(&op4, &arg_2el);

        argument = Transpose::new(arg_2el.clone(), hermitian_transpose).unwrap();
        let op5 = Transpose::new(argument, standard_transpose).unwrap();
        assert_ne!(&op1, &op5);
        assert_eq!(&op5, &Conjugate::new(arg_2el.clone()).unwrap());

        let coef = make_symbol(4u32);
        let arg_wfn = make_wfn_parameter("");
        argument =
            MatrixMul::new(::std::vec![coef.clone(), arg_2el.clone(), arg_wfn.clone()]).unwrap();
        let op6 = Transpose::new(argument, hermitian_transpose).unwrap();
        assert_ne!(&op1, &op6);

        let mat_mul = downcast_from_arc::<MatrixMul>(&op6).unwrap();

        assert_eq!(mat_mul.coefficient(), &Conjugate::new(coef).unwrap());

        assert_eq!(
            mat_mul.factors(),
            vec![
                Transpose::new(
                    MatrixMul::new(::std::vec![arg_2el, arg_wfn]).unwrap(),
                    hermitian_transpose,
                )
                .unwrap()
            ]
        );
    }
}
