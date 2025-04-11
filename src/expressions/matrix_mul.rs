use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Mul, Number, ZeroOperator};
use crate::utils::{
    downcast_from_arc, downcast_from_ref, intern_expr, is_expr_type, is_one_expr, is_zero_expr,
};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct MatrixMul {
    coefficient: Arc<dyn Expr>,
    factors: Vec<Arc<dyn Expr>>,
}

impl MatrixMul {
    // Rules:
    // - Flatten nested MatrixMul
    // - Remove redundant MatrixMul([A]) -> A
    // - Remove empty MatrixMul([]) -> op(0)
    // - Numeric simplifications, e.g. 3 * 2 -> 6, (1/2) * 4 -> 2, 3 * 0 -> 0
    // - Identities, ensure A * 0 = op(0), A * op(0) = op(0)
    // - No polynomial multiplication and expansion, e.g. keeping (A + B) * 2 as-is
    // - Series multiplication, e.g. ((2 * A) * (3 * A)) * A -> 6 * A * A * A
    // - Order of input terms is preserved
    pub fn new(terms: Vec<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        if terms.is_empty() {
            return Ok(ZeroOperator::new());
        }

        let mut all_coefficients = Vec::new();
        let mut all_factors = Vec::new();

        for term in &terms {
            // We may introduce diagonal matrix and dense matrix later, and we
            // may need to check if their dimensions match
            if term.is_scalar() {
                if let Some(num) = downcast_from_arc::<Number>(term) {
                    if num.is_zero() {
                        return Ok(ZeroOperator::new());
                    }
                }
                all_coefficients.push(term.clone());
            } else if is_expr_type::<ZeroOperator>(term) {
                return Ok(ZeroOperator::new());
            } else if let Some(matmul) = downcast_from_arc::<MatrixMul>(term) {
                all_coefficients.push(matmul.coefficient().clone());
                all_factors.extend_from_slice(matmul.factors());
            } else {
                all_factors.push(term.clone());
            }
        }

        if all_coefficients.is_empty() {
            return match all_factors.len() {
                0 => Ok(ZeroOperator::new()), // Should be unreachable
                1 => Ok(all_factors.pop().unwrap()),
                _ => Ok(intern_expr(Arc::new(Self {
                    coefficient: Number::one(),
                    factors: all_factors,
                }))),
            };
        }

        let coefficient = Mul::new(all_coefficients)?;
        match all_factors.len() {
            // Return a pure scalar expression
            0 => Ok(coefficient),
            1 if is_one_expr(&coefficient) => Ok(all_factors.pop().unwrap()),
            _ => Ok(intern_expr(Arc::new(Self {
                coefficient,
                factors: all_factors,
            }))),
        }
    }

    #[inline]
    pub fn coefficient(&self) -> &Arc<dyn Expr> {
        &self.coefficient
    }

    #[inline]
    pub fn factors(&self) -> &[Arc<dyn Expr>] {
        &self.factors
    }
}

impl_mul_traits!(MatrixMul, false);
