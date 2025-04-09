use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::expressions::{Conjugate, HermitianTranspose, MatrixMul, ZeroOperator};
use crate::perturbations::Perturbation;
use crate::utils::{downcast_expr, intern, invalid_expression_error, is_one_expr};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Transpose {
    argument: Arc<dyn Expr>,
}

impl Transpose {
    pub fn new(argument: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if argument.is_scalar() {
            return Err(invalid_expression_error(
                "Transpose::new() - argument must be non-scalar",
                &argument,
            ));
        }

        if argument.is::<ZeroOperator>() {
            return Ok(argument);
        } else if let Some(conj) = downcast_expr::<Conjugate>(&argument) {
            return HermitianTranspose::new(conj.argument().clone());
        } else if let Some(trans) = downcast_expr::<Transpose>(&argument) {
            return Ok(trans.argument().clone());
        } else if let Some(herm) = downcast_expr::<HermitianTranspose>(&argument) {
            return Conjugate::new(herm.argument().clone());
        } else if let Some(matmul) = downcast_expr::<MatrixMul>(&argument) {
            if is_one_expr(matmul.coefficient()) {
                return Ok(intern(Arc::new(Self { argument })));
            }

            let new_arg = MatrixMul::new(matmul.factors().to_vec())?;
            return MatrixMul::new(vec![
                matmul.coefficient().clone(),
                intern(Arc::new(Self { argument: new_arg })),
            ]);
        }

        Ok(intern(Arc::new(Self { argument })))
    }

    #[inline]
    pub fn argument(&self) -> &Arc<dyn Expr> {
        &self.argument
    }
}

impl_unary_expr_traits!(Transpose, false, "{arg}^T");
