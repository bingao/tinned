use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::expressions::{
    Add, Conjugate, HermitianTranspose, MatrixAdd, MatrixMul, Mul, Number, Transpose, ZeroOperator,
};
use crate::perturbations::Perturbation;
use crate::utils::{downcast_expr, intern, invalid_expression_error, is_one_expr};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Trace {
    argument: Arc<dyn Expr>,
}

impl Trace {
    pub fn new(expr: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if expr.is_scalar() {
            return Err(invalid_expression_error("Trace::new()", &expr));
        }

        if expr.is::<ZeroOperator>() {
            Ok(0.into())
        } else if let Some(matrix_add) = downcast_expr::<MatrixAdd>(&expr) {
            let mut terms = Vec::new();
            for term in matrix_add.terms() {
                terms.push(Self::new(term.clone())?);
            }
            Add::new(terms)
        } else if let Some(matrix_mul) = downcast_expr::<MatrixMul>(&expr) {
            let coef = matrix_mul.coefficient();
            let mut factors = matrix_mul.factors().to_vec();

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
            let result = intern(Arc::new(Self { argument: new_mul }));

            if is_one_expr(coef) {
                Ok(result)
            } else {
                Mul::new(vec![coef.clone(), result])
            }
        } else if let Some(conj) = downcast_expr::<Conjugate>(&expr) {
            Conjugate::new(intern(Arc::new(Self { argument: conj.argument().clone() })))
        } else if let Some(trans) = downcast_expr::<Transpose>(&expr) {
            Ok(intern(Arc::new(Self { argument: trans.argument().clone() })))
        } else if let Some(herm) = downcast_expr::<HermitianTranspose>(&expr) {
            Conjugate::new(intern(Arc::new(Self { argument: herm.argument().clone() })))
        } else {
            Ok(intern(Arc::new(Self { argument: expr })))
        }
    }

    #[inline]
    pub fn argument(&self) -> &Arc<dyn Expr> {
        &self.argument
    }
}

impl_unary_expr_traits!(Trace, true, "tr({arg})");
