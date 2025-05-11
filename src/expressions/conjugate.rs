use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{
    Add, DotProduct, HermitianTranspose, MatrixMul, Mul, Number, Power, Transpose, ZeroOperator,
};
use crate::internal::intern_expr;
use crate::perturbations::Perturbation;
use crate::public::{
    NumberTolerance, downcast_from_arc, downcast_from_ref, generic_expression_error, is_expr_type,
    is_one_expr,
};

/// Represents complex conjugation of an expression.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Conjugate {
    argument: Arc<dyn Expr>,
}

impl Conjugate {
    pub fn new(argument: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if let Some(num) = downcast_from_arc::<Number>(&argument) {
            return Ok(num.conjugate().into());
        } else if let Some(add) = downcast_from_arc::<Add>(&argument) {
            let terms: Vec<Arc<dyn Expr>> =
                add.terms().iter().map(|t| Self::new(t.clone())).collect::<Result<_, _>>()?;
            return Add::new(terms);
        } else if let Some(mul) = downcast_from_arc::<Mul>(&argument) {
            let mut new_terms: Vec<Arc<dyn Expr>> =
                mul.factors().iter().map(|f| Self::new(f.clone())).collect::<Result<_, _>>()?;
            new_terms.push(mul.coefficient().conjugate().into());
            return Mul::new(new_terms);
        } else if let Some(power) = downcast_from_arc::<Power>(&argument) {
            let new_base = Self::new(power.base().clone())?;
            return Power::new(new_base, power.exponent());
        } else if is_expr_type::<ZeroOperator>(&argument) {
            return Ok(argument);
        } else if let Some(conj) = downcast_from_arc::<Conjugate>(&argument) {
            return Ok(conj.argument.clone());
        } else if let Some(product) = downcast_from_arc::<DotProduct>(&argument) {
            return product.conjugate();
        } else if let Some(trans) = downcast_from_arc::<Transpose>(&argument) {
            return HermitianTranspose::new(trans.argument().clone());
        } else if let Some(herm) = downcast_from_arc::<HermitianTranspose>(&argument) {
            return Transpose::new(herm.argument().clone());
        } else if let Some(matmul) = downcast_from_arc::<MatrixMul>(&argument) {
            if is_one_expr(matmul.coefficient(), None) {
                return Ok(intern_expr(Arc::new(Self {
                    argument,
                })));
            }

            let new_factors = if matmul.factors().len() == 1 {
                Self::new(Arc::clone(&matmul.factors()[0]))?
            } else {
                intern_expr(Arc::new(Self {
                    argument: MatrixMul::new(matmul.factors().to_vec())?,
                }))
            };
            return MatrixMul::new(vec![Self::new(matmul.coefficient().clone())?, new_factors]);
        }

        Ok(intern_expr(Arc::new(Self {
            argument,
        })))
    }

    #[inline]
    pub fn argument(&self) -> &Arc<dyn Expr> {
        &self.argument
    }
}

impl_unary_expr_traits!(Conjugate, Argument, "conj({arg})");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::number::test_utils::{
        make_number_complex, make_number_f64, make_number_i64, make_number_rational,
    };
    use crate::expressions::symbol::test_utils::make_symbol;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::expressions::{Power, Symbol};
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::is_zero_expr;
    use num_complex::Complex64;
    use num_rational::Rational64;

    test_struct_safety!(Conjugate);

    test_thread_interning!({
        Conjugate::new(
            Mul::new(vec![
                Number::one(),
                Number::from_f64(3.14),
                Number::from_complex(Complex64::new(0.0, -1.0)),
                Number::from_rational(Rational64::new(22, 7)),
                Symbol::new("x"),
                Power::new(Symbol::new("y"), 2).unwrap(),
            ])
            .unwrap(),
        )
        .unwrap()
    });

    #[test]
    fn test_impl_expr() {
        let s1 = make_symbol(4u32);
        let op1 = Conjugate::new(s1.clone()).unwrap();

        let op = downcast_from_arc::<Conjugate>(&op1).unwrap();
        assert_eq!(
            op,
            &Conjugate {
                argument: s1.clone()
            }
        );
        assert_eq!(op.argument(), &s1);

        assert_eq!(op1.hash_key(), format!("Conjugate({})", s1.hash_key()));
        assert_eq!(op1.is_scalar(), s1.is_scalar());
        assert_eq!(format!("{}", op1), format!("conj({})", s1));

        let s2 = make_symbol(8u32);
        let op2 = Conjugate::new(s1.clone()).unwrap();
        let op3 = Conjugate::new(s2.clone()).unwrap();

        assert_eq!(&op1, &op2);
        assert_ne!(&op1, &op3);

        let c1 = make_number_complex(64u32);
        let mut argument = Add::new(vec![c1.clone(), s1.clone(), s2.clone()]).unwrap();
        let conj_add = Conjugate::new(argument.clone()).unwrap();

        assert!(conj_add.is_scalar());
        assert_eq!(
            &conj_add,
            &Add::new(vec![
                Conjugate::new(c1.clone()).unwrap(),
                Conjugate::new(s1.clone()).unwrap(),
                Conjugate::new(s2.clone()).unwrap()
            ])
            .unwrap()
        );
        assert_eq!(&Conjugate::new(conj_add.clone()).unwrap(), &argument);

        argument = Mul::new(vec![c1.clone(), s1.clone(), s2.clone()]).unwrap();
        let conj_mul = Conjugate::new(argument.clone()).unwrap();

        assert!(conj_mul.is_scalar());
        assert_eq!(
            &conj_mul,
            &Mul::new(vec![
                Conjugate::new(c1.clone()).unwrap(),
                Conjugate::new(s1.clone()).unwrap(),
                Conjugate::new(s2.clone()).unwrap()
            ])
            .unwrap()
        );
        assert_eq!(&Conjugate::new(conj_mul.clone()).unwrap(), &argument);

        let exponent: i64 = rand::random_range(2..=16);
        argument = Power::new(s1.clone(), exponent).unwrap();
        let conj_power = Conjugate::new(argument.clone()).unwrap();

        assert!(conj_power.is_scalar());
        assert_eq!(
            &conj_power,
            &Power::new(Conjugate::new(s1.clone()).unwrap(), exponent).unwrap()
        );
        assert_eq!(&Conjugate::new(conj_power.clone()).unwrap(), &argument);

        let op_a = make_wfn_parameter("");
        let op_b = make_wfn_parameter("");
        argument =
            MatrixMul::new(vec![c1.clone(), s1.clone(), op_a.clone(), op_b.clone()]).unwrap();
        let conj_matmul = Conjugate::new(argument.clone()).unwrap();

        assert!(!conj_matmul.is_scalar());
        assert_eq!(
            &conj_matmul,
            &MatrixMul::new(vec![
                Conjugate::new(Mul::new(vec![c1.clone(), s1.clone()]).unwrap()).unwrap(),
                Conjugate::new(MatrixMul::new(vec![op_a.clone(), op_b.clone()]).unwrap()).unwrap()
            ])
            .unwrap()
        );
        assert_eq!(&Conjugate::new(conj_matmul.clone()).unwrap(), &argument);

        argument = DotProduct::new(op_a.clone(), true, op_b.clone(), false).unwrap();
        let mut conj_dot = Conjugate::new(argument.clone()).unwrap();

        assert!(conj_dot.is_scalar());
        assert_eq!(
            &conj_dot,
            &DotProduct::new(
                Conjugate::new(op_a.clone()).unwrap(),
                true,
                Conjugate::new(op_b.clone()).unwrap(),
                false
            )
            .unwrap()
        );
        assert_eq!(&Conjugate::new(conj_dot.clone()).unwrap(), &argument);

        argument = DotProduct::new(op_a.clone(), false, op_b.clone(), false).unwrap();
        conj_dot = Conjugate::new(argument.clone()).unwrap();

        assert!(conj_dot.is_scalar());
        assert_eq!(
            &conj_dot,
            &DotProduct::new(
                Conjugate::new(op_a.clone()).unwrap(),
                false,
                Conjugate::new(op_b.clone()).unwrap(),
                false
            )
            .unwrap()
        );
        assert_eq!(&Conjugate::new(conj_dot.clone()).unwrap(), &argument);

        argument = make_wfn_parameter("");
        let conj_dagger =
            Conjugate::new(HermitianTranspose::new(argument.clone()).unwrap()).unwrap();
        let conj_trans = Conjugate::new(Transpose::new(argument.clone()).unwrap()).unwrap();

        assert!(!conj_dagger.is_scalar());
        assert!(!conj_trans.is_scalar());
        assert_eq!(&conj_dagger, &Transpose::new(argument.clone()).unwrap());
        assert_eq!(&conj_trans, &HermitianTranspose::new(argument.clone()).unwrap());
        assert_eq!(
            &Conjugate::new(conj_dagger.clone()).unwrap(),
            &HermitianTranspose::new(argument.clone()).unwrap()
        );
        assert_eq!(
            &Conjugate::new(conj_trans.clone()).unwrap(),
            &Transpose::new(argument.clone()).unwrap()
        );

        assert_eq!(&Conjugate::new(ZeroOperator::new()).unwrap(), &ZeroOperator::new());
    }

    #[test]
    fn test_differentiation() {
        let mut argument = make_wfn_parameter("");
        let mut op = Conjugate::new(argument.clone()).unwrap();
        let p = make_perturbation_symbol(4u32, 4u32);

        assert_eq!(
            &op.differentiate(&p).unwrap(),
            &Conjugate::new(argument.differentiate(&p).unwrap()).unwrap()
        );

        argument =
            DotProduct::new(make_wfn_parameter(""), true, make_wfn_parameter(""), false).unwrap();
        op = Conjugate::new(argument.clone()).unwrap();

        assert_eq!(
            &op.differentiate(&p).unwrap(),
            &Conjugate::new(argument.differentiate(&p).unwrap()).unwrap()
        );

        assert_eq!(
            &Conjugate::new(make_symbol(4u32)).unwrap().differentiate(&p).unwrap(),
            &Number::zero()
        );
    }

    #[test]
    fn test_serialization() {
        let op = Conjugate::new(
            Mul::new(vec![
                make_number_i64(256u32),
                make_number_f64(64u32),
                make_number_complex(64u32),
                make_number_rational(256u32),
                make_symbol(4u32),
                Power::new(make_symbol(4u32), rand::random_range(-256..=256) as i64).unwrap(),
            ])
            .unwrap(),
        )
        .unwrap();
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let s1 = make_symbol(4u32);
        let op1 = Conjugate::new(s1.clone()).unwrap();

        assert!(is_expr_type::<Conjugate>(&op1));
        assert!(!is_zero_expr(&op1, None));
        assert!(!is_one_expr(&op1, None));

        let s2 = make_symbol(8u32);
        let op2 = Conjugate::new(s1.clone()).unwrap();
        let op3 = Conjugate::new(s2.clone()).unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));
    }
}
