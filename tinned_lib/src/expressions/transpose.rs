use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{Conjugate, HermitianTranspose, MatrixMul, ZeroOperator};
use crate::internal::intern_expr;
use crate::perturbations::Perturbation;
use crate::public::{
    NumberTolerance, downcast_from_arc, downcast_from_ref, expression_error,
    generic_expression_error, is_expr_type, is_one_expr,
};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Transpose {
    argument: Arc<dyn Expr>,
}

impl Transpose {
    pub fn new(argument: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
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
            return HermitianTranspose::new(conj.argument().clone());
        } else if let Some(trans) = downcast_from_arc::<Transpose>(&argument) {
            return Ok(trans.argument().clone());
        } else if let Some(herm) = downcast_from_arc::<HermitianTranspose>(&argument) {
            return Conjugate::new(herm.argument().clone());
        } else if let Some(mat_mul) = downcast_from_arc::<MatrixMul>(&argument) {
            if is_one_expr(mat_mul.coefficient(), None) {
                return Ok(intern_expr(Arc::new(Self {
                    argument,
                })));
            }

            let new_arg = MatrixMul::new(mat_mul.factors().to_vec())?;
            return MatrixMul::new(vec![
                mat_mul.coefficient().clone(),
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

impl_unary_expr_traits!(Transpose, False, "{arg}^T");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::symbol::test_utils::make_symbol;
    use crate::expressions::two_elec_operator::test_utils::make_two_elec_operator;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;

    test_transpose!(Transpose, HermitianTranspose, false, "{arg}^T");
}
