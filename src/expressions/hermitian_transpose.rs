use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Conjugate, MatrixMul, Transpose, ZeroOperator};
use crate::utils::{
    downcast_from_arc, downcast_from_ref, expression_error, generic_expression_error, intern_expr,
    is_expr_type, is_one_expr,
};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct HermitianTranspose {
    argument: Arc<dyn Expr>,
}

impl HermitianTranspose {
    pub fn new(argument: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if argument.is_scalar() {
            return Err(expression_error(
                "HermitianTranspose::new() - argument must be non-scalar",
                &argument,
                None,
            ));
        }

        if is_expr_type::<ZeroOperator>(&argument) {
            return Ok(argument);
        } else if let Some(conj) = downcast_from_arc::<Conjugate>(&argument) {
            return Transpose::new(conj.argument().clone());
        } else if let Some(trans) = downcast_from_arc::<Transpose>(&argument) {
            return Conjugate::new(trans.argument().clone());
        } else if let Some(herm) = downcast_from_arc::<HermitianTranspose>(&argument) {
            return Ok(herm.argument().clone());
        } else if let Some(matmul) = downcast_from_arc::<MatrixMul>(&argument) {
            if is_one_expr(matmul.coefficient(), None) {
                return Ok(intern_expr(Arc::new(Self {
                    argument,
                })));
            }

            let new_arg = MatrixMul::new(matmul.factors().to_vec())?;
            return MatrixMul::new(vec![
                Conjugate::new(matmul.coefficient().clone())?,
                intern_expr(Arc::new(Self {
                    argument: new_arg,
                })),
            ]);
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

impl_unary_expr_traits!(HermitianTranspose, false, "{arg}^H");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::symbol::test_utils::make_symbol;
    use crate::expressions::two_elec_operator::test_utils::make_two_elec_operator;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;

    test_transpose!(HermitianTranspose, Transpose, true, "{arg}^H");
}
