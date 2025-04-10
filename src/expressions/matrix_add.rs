use std::collections::HashMap;
use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Add, MatrixMul, Number, ZeroOperator};
use crate::utils::{
    downcast_from_arc, downcast_from_ref, intern_expr, invalid_expression_error, is_expr_type,
    is_one_expr, is_zero_expr, unreachable_error,
};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct MatrixAdd {
    terms: Vec<Arc<dyn Expr>>,
}

impl MatrixAdd {
    // Rules:
    // - Flatten nested MatrixAdd: (A + B) + (A + C) -> 2A + B + C
    // - Remove redundant MatrixAdd([A]) -> A
    // - Remove empty MatrixAdd([]) -> op(0)
    // - Combine like terms: 2*A*B + 3*A*B -> 5*A*B
    // - Identities: A + op(0) = A
    // - Sort terms based on hash values
    pub fn new(terms: Vec<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        let mut merged: HashMap<u64, (Arc<dyn Expr>, Vec<Arc<dyn Expr>>)> = HashMap::new();

        #[inline]
        fn collect_terms(
            expr: &Arc<dyn Expr>,
            merged: &mut HashMap<u64, (Arc<dyn Expr>, Vec<Arc<dyn Expr>>)>,
        ) -> Result<(), TinnedError> {
            if expr.is_scalar() {
                return Err(invalid_expression_error("MatrixAdd::new()", &expr));
            }

            if is_expr_type::<ZeroOperator>(expr) {
                return Ok(());
            } else if let Some(matadd) = downcast_from_arc::<MatrixAdd>(expr) {
                for term in matadd.terms() {
                    collect_terms(term, merged)?;
                }
            } else if let Some(matmul) = downcast_from_arc::<MatrixMul>(expr) {
                if matmul.factors().is_empty() {
                    return Err(unreachable_error(
                        "MatrixAdd::new() got MatrixMul with empty factors",
                        &expr,
                    ));
                }
                let base_expr = if matmul.factors().len() == 1 {
                    matmul.factors()[0].clone()
                } else {
                    MatrixMul::new(matmul.factors().to_vec())?
                };
                let key = base_expr.fast_hash();
                let coef = matmul.coefficient();
                if let Some((existing_expr, existing_coef)) = merged.get_mut(&key) {
                    if existing_expr == &base_expr {
                        existing_coef.push(coef.clone());
                        return Ok(());
                    }
                }
                merged.insert(key, (base_expr, vec![coef.clone()]));
            } else {
                let key = expr.fast_hash();
                if let Some((existing_expr, existing_coef)) = merged.get_mut(&key) {
                    if existing_expr == expr {
                        existing_coef.push(Arc::new(Number::Integer(1)));
                        return Ok(());
                    }
                }
                merged.insert(key, (expr.clone(), vec![Arc::new(Number::Integer(1))]));
            }

            Ok(())
        }

        for term in &terms {
            collect_terms(term, &mut merged)?;
        }

        let mut simplified_terms = Vec::new();

        for (expr, all_coef) in merged.into_values() {
            let coef = Add::new(all_coef)?;
            if !is_zero_expr(&coef) {
                if is_one_expr(&coef) {
                    simplified_terms.push(expr);
                } else {
                    simplified_terms.push(MatrixMul::new(vec![coef, expr])?);
                }
            }
        }

        match simplified_terms.len() {
            0 => Ok(ZeroOperator::new()),
            1 => Ok(simplified_terms.pop().unwrap()),
            _ => {
                simplified_terms.sort_by_key(|term| term.fast_hash());
                Ok(intern_expr(Arc::new(Self { terms: simplified_terms })))
            },
        }
    }

    #[inline]
    pub fn terms(&self) -> &[Arc<dyn Expr>] {
        &self.terms
    }
}

impl_add_traits!(MatrixAdd, false);
