use std::collections::HashMap;
use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Mul, Number};
use crate::utils::{
    downcast_from_arc, downcast_from_ref, intern_expr, invalid_expression_error, unreachable_error,
};

// Addition Expression
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Add {
    terms: Vec<Arc<dyn Expr>>,
}

impl Add {
    // Rules:
    // - Flatten nested Add: (x + 2) + (x + 3) -> 2x + 5
    // - Remove redundant Add([x]) -> x
    // - Remove empty Add([]) -> 0
    // - Numeric simplifications: 3 + 5 -> 8
    // - Combine like terms: 2*x*y + 3*x*y -> 5*x*y
    // - Identities: x + 0 = x
    // - Sort terms based on hash values
    pub fn new(terms: Vec<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        let mut constant = Number::Integer(0);
        // Key: fast hash of term (u64), Value: (expr, accumulated coefficient)
        let mut merged: HashMap<u64, (Arc<dyn Expr>, Number)> = HashMap::new();

        #[inline]
        fn collect_terms(
            expr: &Arc<dyn Expr>,
            constant: &mut Number,
            merged: &mut HashMap<u64, (Arc<dyn Expr>, Number)>,
        ) -> Result<(), TinnedError> {
            if !expr.is_scalar() {
                return Err(invalid_expression_error("Add::new()", &expr));
            }

            if let Some(num) = downcast_from_arc::<Number>(expr) {
                *constant = constant.add(num);
            } else if let Some(add) = downcast_from_arc::<Add>(expr) {
                for term in add.terms() {
                    collect_terms(term, constant, merged)?;
                }
            } else if let Some(mul) = downcast_from_arc::<Mul>(expr) {
                if mul.factors().is_empty() {
                    return Err(unreachable_error("Add::new() got Mul with empty factors", expr));
                }

                let base_expr = if mul.factors().len() == 1 {
                    mul.factors()[0].clone()
                } else {
                    Mul::new(mul.factors().to_vec())?
                };
                let key = base_expr.fast_hash();
                let coef = mul.coefficient();

                if let Some((existing_expr, existing_coef)) = merged.get_mut(&key) {
                    if existing_expr == &base_expr {
                        *existing_coef = existing_coef.add(coef);
                        return Ok(());
                    }
                }
                merged.insert(key, (base_expr, coef.clone()));
            } else {
                let key = expr.fast_hash();
                if let Some((existing_expr, existing_coef)) = merged.get_mut(&key) {
                    if existing_expr == expr {
                        *existing_coef = existing_coef.add(&Number::Integer(1));
                        return Ok(());
                    }
                }
                merged.insert(key, (expr.clone(), Number::Integer(1)));
            }

            Ok(())
        }

        for term in &terms {
            collect_terms(term, &mut constant, &mut merged)?;
        }

        let mut simplified_terms: Vec<Arc<dyn Expr>> = Vec::with_capacity(merged.len());

        for (expr, coef) in merged.into_values() {
            if !coef.is_zero() {
                if coef.is_one() {
                    simplified_terms.push(expr);
                } else {
                    simplified_terms.push(Mul::new(vec![coef.into(), expr])?);
                }
            }
        }

        if !constant.is_zero() {
            simplified_terms.push(intern_expr(Arc::new(constant)));
        }

        match simplified_terms.len() {
            0 => Ok(Number::zero()),
            1 => Ok(simplified_terms.pop().unwrap()),
            _ => {
                simplified_terms.sort_by_key(|term| term.fast_hash());
                Ok(intern_expr(Arc::new(Self {
                    terms: simplified_terms,
                })))
            },
        }
    }

    #[inline]
    pub fn terms(&self) -> &[Arc<dyn Expr>] {
        &self.terms
    }
}

impl_add_traits!(Add, true);
