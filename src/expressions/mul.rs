use std::collections::HashMap;
use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Number, Power};
use crate::utils::{
    downcast_from_arc, downcast_from_ref, intern, invalid_expression_error, is_zero_expr,
};

// Multiplication Expression
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Mul {
    coefficient: Number,
    factors: Vec<Arc<dyn Expr>>,
}

impl Mul {
    // Rules:
    // - Flatten nested Mul
    // - Remove redundant Mul([x]) -> x
    // - Remove empty Mul([]) -> 0
    // - Numeric simplifications, e.g. 3 * 2 -> 6, (1/2) * 4 -> 2, 3 * 0 -> 0
    // - Merge powers in Mul, e.g. x^a * x^b -> x^(a+b)
    // - Ensure x * x^a = x^(a+1), x^a * x = x^(a+1)
    // - Identities, ensure x * 1 = x, x * 0 = 0
    // - No polynomial multiplication and expansion, e.g. keeping (x + y) * 2 as-is
    // - Series multiplication, e.g. ((2 * x) * (3 * x)) * x -> 6 * x^3
    // - Sort terms based on hash values
    pub fn new(terms: Vec<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        if terms.is_empty() {
            return Ok(0.into());
        }

        let mut coefficient = Number::Integer(1);
        let mut power_map: HashMap<u64, (Arc<dyn Expr>, i64)> = HashMap::new();

        #[inline]
        fn collect_terms(
            expr: &Arc<dyn Expr>,
            coefficient: &mut Number,
            power_map: &mut HashMap<u64, (Arc<dyn Expr>, i64)>,
        ) -> Result<bool, TinnedError> {
            if !expr.is_scalar() {
                return Err(invalid_expression_error("Mul::new()", &expr));
            }

            if let Some(mul) = downcast_from_arc::<Mul>(expr) {
                *coefficient = coefficient.mul(mul.coefficient());
                for factor in mul.factors() {
                    if collect_terms(factor, coefficient, power_map)? {
                        return Ok(true);
                    }
                }
            } else if let Some(num) = downcast_from_arc::<Number>(expr) {
                if num.is_zero() {
                    return Ok(true); // Multiplication by 0 -> entire result is zero
                } else if !num.is_one() {
                    *coefficient = coefficient.mul(num);
                }
            } else if let Some(pow) = downcast_from_arc::<Power>(expr) {
                let key = pow.base().fast_hash();
                let base = pow.base().clone();
                let exp = pow.exponent();
                let entry = power_map.entry(key).or_insert((base, 0));
                entry.1 += exp;
            } else {
                let key = expr.fast_hash();
                let entry = power_map.entry(key).or_insert((expr.clone(), 0));
                entry.1 += 1;
            }

            Ok(false)
        }

        for term in &terms {
            if collect_terms(term, &mut coefficient, &mut power_map)? {
                return Ok(0.into());
            }
        }

        let mut simplified_factors: Vec<Arc<dyn Expr>> = Vec::with_capacity(power_map.len());

        for (base, exp) in power_map.into_values() {
            if exp == 0 {
                continue;
            } else if exp == 1 {
                simplified_factors.push(base);
            } else {
                simplified_factors.push(Power::new(base, exp)?);
            }
        }

        simplified_factors.sort_by_key(|f| f.fast_hash());

        match simplified_factors.len() {
            0 => Ok(coefficient.into()),
            1 if coefficient.is_one() => Ok(simplified_factors.pop().unwrap()),
            _ => Ok(intern(Arc::new(Self { coefficient, factors: simplified_factors }))),
        }
    }

    #[inline]
    pub fn coefficient(&self) -> &Number {
        &self.coefficient
    }

    #[inline]
    pub fn factors(&self) -> &[Arc<dyn Expr>] {
        &self.factors
    }
}

impl_mul_traits!(Mul, true);
