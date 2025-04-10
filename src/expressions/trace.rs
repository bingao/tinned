use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{
    Add, Conjugate, HermitianTranspose, MatrixAdd, MatrixMul, Mul, Number, Transpose, ZeroOperator,
};
use crate::utils::{
    downcast_from_arc, downcast_from_ref, intern_expr, invalid_expression_error, is_expr_type,
    is_one_expr,
};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Trace {
    argument: Arc<dyn Expr>,
}

impl Trace {
    pub fn new(argument: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if argument.is_scalar() {
            return Err(invalid_expression_error("Trace::new()", &argument));
        }

        if is_expr_type::<ZeroOperator>(&argument) {
            Ok(Number::zero())
        } else if let Some(matadd) = downcast_from_arc::<MatrixAdd>(&argument) {
            let mut terms = Vec::new();
            for term in matadd.terms() {
                terms.push(Self::new(term.clone())?);
            }
            Add::new(terms)
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
            let result = intern_expr(Arc::new(Self { argument: new_mul }));

            if is_one_expr(coef) {
                Ok(result)
            } else {
                Mul::new(vec![coef.clone(), result])
            }
        } else if let Some(conj) = downcast_from_arc::<Conjugate>(&argument) {
            Conjugate::new(intern_expr(Arc::new(Self { argument: conj.argument().clone() })))
        } else if let Some(trans) = downcast_from_arc::<Transpose>(&argument) {
            Ok(intern_expr(Arc::new(Self { argument: trans.argument().clone() })))
        } else if let Some(herm) = downcast_from_arc::<HermitianTranspose>(&argument) {
            Conjugate::new(intern_expr(Arc::new(Self { argument: herm.argument().clone() })))
        } else {
            Ok(intern_expr(Arc::new(Self { argument })))
        }
    }

    #[inline]
    pub fn argument(&self) -> &Arc<dyn Expr> {
        &self.argument
    }
}

impl_unary_expr_traits!(Trace, true, "tr({arg})");
