use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{
    Add, Conjugate, HermitianTranspose, MatrixAdd, MatrixMul, Mul, Number, Transpose, ZeroOperator,
};
use crate::internal::intern_expr;
use crate::perturbations::Perturbation;
use crate::public::{
    NumberTolerance, downcast_from_arc, downcast_from_ref, expression_error,
    generic_expression_error, is_expr_type, is_one_expr,
};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Trace {
    argument: Arc<dyn Expr>,
}

impl Trace {
    pub fn new(argument: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if argument.is_scalar() {
            return Err(expression_error("Trace::new() got a scalar argument", &argument, None));
        }

        if is_expr_type::<ZeroOperator>(&argument) {
            Ok(Number::zero())
        } else if let Some(matadd) = downcast_from_arc::<MatrixAdd>(&argument) {
            let terms = matadd.terms();
            let mut new_terms = Vec::with_capacity(terms.len());
            for term in terms {
                new_terms.push(Self::new(term.clone())?);
            }
            Add::new(new_terms)
        } else if let Some(matmul) = downcast_from_arc::<MatrixMul>(&argument) {
            let coef = matmul.coefficient();
            let mut factors = matmul.factors().to_vec();

            if factors.len() > 1 {
                // circular shift to bring minimal hash to front
                let min_idx = factors
                    .iter()
                    .enumerate()
                    .min_by_key(|(_, f)| f.fast_hash())
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                factors.rotate_left(min_idx);
            }

            let new_mul = MatrixMul::new(factors)?;
            let result = intern_expr(Arc::new(Self {
                argument: new_mul,
            }));

            if is_one_expr(coef, None) {
                Ok(result)
            } else {
                Mul::new(vec![coef.clone(), result])
            }
        } else if let Some(conj) = downcast_from_arc::<Conjugate>(&argument) {
            Conjugate::new(intern_expr(Arc::new(Self {
                argument: conj.argument().clone(),
            })))
        } else if let Some(trans) = downcast_from_arc::<Transpose>(&argument) {
            Ok(intern_expr(Arc::new(Self {
                argument: trans.argument().clone(),
            })))
        } else if let Some(herm) = downcast_from_arc::<HermitianTranspose>(&argument) {
            Conjugate::new(intern_expr(Arc::new(Self {
                argument: herm.argument().clone(),
            })))
        } else {
            Ok(intern_expr(Arc::new(Self {
                argument,
            })))
        }
    }

    #[inline]
    pub fn argument(&self) -> &Arc<dyn Expr> {
        &self.argument
    }
}

impl_unary_expr_traits!(Trace, True, "tr({arg})");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::number::test_utils::make_number_complex;
    use crate::expressions::two_elec_operator::test_utils::make_two_elec_operator;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;

    test_unary_oper_properties!(Trace);

    #[test]
    fn test_impl_expr() {
        let op0 = Trace::new(ZeroOperator::new()).unwrap();
        assert!(crate::public::is_zero_expr(&op0, None));

        let arg_2el = make_two_elec_operator("", None);
        let op1 = Trace::new(arg_2el.clone()).unwrap();

        let op = downcast_from_arc::<Trace>(&op1).unwrap();
        assert_eq!(
            op,
            &Trace {
                argument: arg_2el.clone()
            }
        );
        assert_eq!(op.argument(), &arg_2el);

        let op2 = Trace::new(arg_2el.clone()).unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        assert_eq!(op1.hash_key(), format!("Trace({})", arg_2el.hash_key()));
        assert!(op1.is_scalar());
        assert_eq!(format!("{}", op1), format!("tr({})", arg_2el));

        let mut argument = Conjugate::new(arg_2el.clone()).unwrap();
        let op3 = Trace::new(argument).unwrap();
        assert_ne!(&op1, &op3);
        assert_eq!(&op3, &Conjugate::new(Trace::new(arg_2el.clone()).unwrap()).unwrap());

        argument = Transpose::new(arg_2el.clone()).unwrap();
        let op4 = Trace::new(argument).unwrap();
        assert!(Arc::ptr_eq(&op1, &op4));
        assert_eq!(&op1, &op4);

        argument = HermitianTranspose::new(arg_2el.clone()).unwrap();
        let op5 = Trace::new(argument).unwrap();
        assert_ne!(&op1, &op5);
        assert_eq!(&op5, &Conjugate::new(Trace::new(arg_2el.clone()).unwrap()).unwrap());

        let arg_wfn = make_wfn_parameter("");
        argument = MatrixAdd::new(vec![arg_2el.clone(), arg_wfn.clone()]).unwrap();
        let op6 = Trace::new(argument).unwrap();
        assert_ne!(&op1, &op6);
        assert_eq!(
            &op6,
            &Add::new(vec![
                Trace::new(arg_2el.clone()).unwrap(),
                Trace::new(arg_wfn.clone()).unwrap()
            ])
            .unwrap()
        );

        let coef = make_number_complex(100u32);
        argument = MatrixMul::new(vec![coef.clone(), arg_2el.clone(), arg_wfn.clone()]).unwrap();
        let op7 = Trace::new(argument).unwrap();
        assert_ne!(&op1, &op7);

        let mul = downcast_from_arc::<Mul>(&op7).unwrap();
        //let expected: Arc<dyn Expr> = mul.coefficient().clone().into();
        let expected: Arc<dyn Expr> = mul.coefficient().into();
        assert_eq!(&expected, &coef);
        let terms = if arg_2el.hash_key() <= arg_wfn.hash_key() {
            vec![arg_2el, arg_wfn]
        } else {
            vec![arg_wfn, arg_2el]
        };
        assert_eq!(mul.factors(), vec![Trace::new(MatrixMul::new(terms).unwrap()).unwrap()]);
    }
}
