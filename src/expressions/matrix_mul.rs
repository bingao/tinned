use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::expressions::{MatrixAdd, Mul, Number, ZeroOperator};
use crate::perturbations::Perturbation;
use crate::utils::{downcast_expr, intern, is_one_expr, is_zero_expr};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
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
                if let Some(num) = downcast_expr::<Number>(term) {
                    if num.is_zero() {
                        return Ok(ZeroOperator::new());
                    }
                }
                all_coefficients.push(term.clone());
            } else if term.is::<ZeroOperator>() {
                return Ok(ZeroOperator::new());
            } else if let Some(mul) = downcast_expr::<MatrixMul>(term) {
                all_coefficients.push(mul.coefficient().clone());
                all_factors.extend_from_slice(mul.factors());
            } else {
                all_factors.push(term.clone());
            }
        }

        if all_coefficients.is_empty() {
            return match all_factors.len() {
                0 => Ok(ZeroOperator::new()), // Should be unreachable
                1 => Ok(all_factors.pop().unwrap()),
                _ => Ok(intern(Arc::new(Self { coefficient: 1.into(), factors: all_factors }))),
            };
        }

        let coefficient = Mul::new(all_coefficients)?;
        match all_factors.len() {
            // Return a pure scalar expression
            0 => Ok(coefficient),
            1 if is_one_expr(&coefficient) => Ok(all_factors.pop().unwrap()),
            _ => Ok(intern(Arc::new(Self { coefficient, factors: all_factors }))),
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

impl Expr for MatrixMul {
    #[inline]
    fn as_any(&self) -> &dyn Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        let keys: Vec<String> = self.factors.iter().map(|f| f.hash_key()).collect();
        format!("MatrixMul({}; {})", self.coefficient.hash_key(), keys.join(","))
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        false
    }

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        if let Some(mul) = downcast_expr::<MatrixMul>(other) {
            self.coefficient == mul.coefficient && self.factors == mul.factors
        } else {
            false
        }
    }

    fn differentiate(&self, s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_coef = self.coefficient.differentiate(s)?;
        let diff_factors: Vec<_> = self.factors.iter().map(|f| f.differentiate(s)?).collect();

        let mut results = Vec::new();

        // If coefficient's derivative is non-zero, append it as one result
        if !is_zero_expr(&diff_coef) {
            let mut new_terms = self.factors.clone();
            new_terms.push(diff_coef);
            results.push(Self::new(new_terms)?);
        }

        for (i, diff) in diff_factors.iter().enumerate() {
            if is_zero_expr(diff) {
                continue;
            }

            let mut new_terms = self.factors.clone();
            new_terms[i] = diff.clone();
            new_terms.push(self.coefficient.clone());

            results.push(Self::new(new_terms)?);
        }

        MatrixAdd::new(results)
    }
}

impl Display for MatrixMul {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        let mut parts = Vec::new();
        if !is_one_expr(&self.coefficient) || self.factors.is_empty() {
            parts.push(format!("{}", self.coefficient));
        }
        for factor in &self.factors {
            parts.push(format!("{}", factor));
        }
        write!(f, "{}", parts.join(" * "))
    }
}
