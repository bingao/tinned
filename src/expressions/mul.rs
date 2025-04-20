use std::collections::HashMap;
use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Number, Power};
use crate::utils::{
    downcast_from_arc, downcast_from_ref, intern_expr, invalid_expression_error, is_zero_expr,
    join_exprs_for_display, join_exprs_for_hash,
};

// Multiplication Expression
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
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
            return Ok(Number::zero());
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
                *coefficient = coefficient.mul(&mul.coefficient);
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
                return Ok(Number::zero());
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
            _ => Ok(intern_expr(Arc::new(Self {
                coefficient,
                factors: simplified_factors,
            }))),
        }
    }

    /// Returns the coefficient as a reference to the internal `Number`.
    #[inline]
    pub fn coefficient(&self) -> &Number {
        &self.coefficient
    }

    #[inline]
    pub fn factors(&self) -> &[Arc<dyn Expr>] {
        &self.factors
    }
}

const DEFAULT_HASH_DELIMITER: &str = ";";
const DEFAULT_FMT_DELIMITER: &str = " * ";

impl_mul_traits!(Mul, DEFAULT_HASH_DELIMITER, DEFAULT_FMT_DELIMITER, true);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::number::test_utils::{
        make_number_complex, make_number_f64, make_number_i64, make_number_rational,
    };
    use crate::expressions::symbol::test_utils::make_symbol;
    use crate::expressions::{Add, Symbol};
    use crate::utils::{is_expr_type, is_one_expr};
    use num_complex::Complex64;
    use num_rational::Rational64;

    test_struct_safety!(Mul);

    test_thread_interning!({
        Mul::new(vec![
            Number::one(),
            Number::from_f64(3.14),
            Number::from_complex(Complex64::new(0.0, -1.0)),
            Number::from_rational(Rational64::new(22, 7)),
            Symbol::new("x"),
            Power::new(Symbol::new("y"), 2).unwrap(),
        ])
        .unwrap()
    });

    #[test]
    fn test_impl_expr() {
        let c1 = make_number_complex(64u32);
        let x = Symbol::new("x");
        let y = Symbol::new("y");
        let z = Symbol::new("z");
        let mul1 = Mul::new(vec![c1.clone(), x.clone(), y.clone(), z.clone()]).unwrap();

        assert!(is_expr_type::<Mul>(&mul1));

        let c1_cast = downcast_from_arc::<Number>(&c1).unwrap();
        let mut mul = downcast_from_arc::<Mul>(&mul1).unwrap();
        let mut asc_factors = vec![x.clone(), y.clone(), z.clone()];
        let mut desc_factors = vec![x.clone(), y.clone(), z.clone()];

        asc_factors.sort_by_key(|f| f.fast_hash());
        desc_factors.sort_by_key(|f| std::cmp::Reverse(f.fast_hash()));

        // - Sort terms based on hash values
        assert_eq!(
            mul,
            &Mul {
                coefficient: c1_cast.clone(),
                factors: asc_factors.clone(),
            }
        );
        assert_eq!(mul.coefficient(), c1_cast);
        assert_eq!(mul.factors(), &asc_factors);
        assert_ne!(mul.factors(), &desc_factors);

        assert_eq!(
            mul1.hash_key(),
            format!(
                "Mul({}{}{})",
                c1_cast.hash_key(),
                DEFAULT_HASH_DELIMITER,
                join_exprs_for_hash(&asc_factors, DEFAULT_HASH_DELIMITER),
            )
        );
        assert!(mul1.is_scalar());
        if c1_cast.is_one() {
            assert_eq!(
                format!("{}", mul1),
                format!("{}", join_exprs_for_display(&asc_factors, DEFAULT_FMT_DELIMITER))
            );
        } else {
            assert_eq!(
                format!("{}", mul1),
                format!(
                    "{}{}{}",
                    c1_cast,
                    DEFAULT_FMT_DELIMITER,
                    join_exprs_for_display(&asc_factors, DEFAULT_FMT_DELIMITER),
                )
            );
        }

        let mul2 = Mul::new(vec![c1.clone(), x.clone(), y.clone(), z.clone()]).unwrap();
        let mul3 = Mul::new(vec![x.clone(), y.clone(), z.clone()]).unwrap();
        let mul4 = Mul::new(vec![c1.clone(), x.clone(), y.clone()]).unwrap();

        assert_eq!(&mul1, &mul2);
        assert_ne!(&mul1, &mul3);
        assert_ne!(&mul1, &mul4);

        mul = downcast_from_arc::<Mul>(&mul3).unwrap();

        assert_eq!(mul.coefficient(), &Number::Integer(1));

        // - Remove empty Mul([]) -> 0
        assert_eq!(&Mul::new(vec![]).unwrap(), &Number::zero());

        // - Remove redundant Mul([x]) -> x
        assert_eq!(&Mul::new(vec![x.clone()]).unwrap(), &x);

        // - Merge powers in Mul, e.g. x^a * x^b -> x^(a+b)
        let exponent1: i64 = rand::random_range(1..=16);
        let exponent2: i64 = rand::random_range(-64..=64);
        assert_eq!(
            &Mul::new(vec![
                Power::new(x.clone(), exponent1).unwrap(),
                Power::new(x.clone(), exponent2).unwrap(),
            ])
            .unwrap(),
            &Power::new(x.clone(), exponent1 + exponent2).unwrap()
        );

        // - Ensure x * x^a = x^(a+1), x^a * x = x^(a+1)
        assert_eq!(
            &Mul::new(vec![x.clone(), Power::new(x.clone(), exponent1).unwrap()]).unwrap(),
            &Power::new(x.clone(), exponent1 + 1).unwrap()
        );
        assert_eq!(
            &Mul::new(vec![Power::new(x.clone(), exponent1).unwrap(), x.clone()]).unwrap(),
            &Power::new(x.clone(), exponent1 + 1).unwrap()
        );

        // - Identities, ensure x * 1 = x, x * 0 = 0
        assert_eq!(&Mul::new(vec![x.clone(), Number::one()]).unwrap(), &x);
        assert_eq!(&Mul::new(vec![x.clone(), Number::zero()]).unwrap(), &Number::zero());

        // - Numeric simplifications, e.g. 3 * 2 -> 6, (1/2) * 4 -> 2, 3 * 0 -> 0
        let c2 = make_number_complex(64u32);
        let c3 = make_number_complex(64u32);
        let mul5 = Mul::new(vec![c1.clone(), c2.clone(), c3.clone()]).unwrap();

        assert!(is_expr_type::<Number>(&mul5));

        let num = downcast_from_arc::<Number>(&mul5).unwrap();
        let c2_cast = downcast_from_arc::<Number>(&c2).unwrap();
        let c3_cast = downcast_from_arc::<Number>(&c3).unwrap();

        assert_eq!(num, &c1_cast.mul(&c2_cast.mul(&c3_cast)));

        // - Flatten nested Mul
        // - Series multiplication, e.g. ((2 * x) * (3 * x)) * x -> 6 * x^3
        assert_eq!(
            &Mul::new(vec![
                Mul::new(vec![
                    Mul::new(vec![c1.clone(), x.clone()]).unwrap(),
                    Mul::new(vec![c2.clone(), x.clone()]).unwrap(),
                ])
                .unwrap(),
                c3.clone(),
                x.clone()
            ])
            .unwrap(),
            &Mul::new(vec![c1.clone(), c2.clone(), c3.clone(), Power::new(x.clone(), 3).unwrap(),])
                .unwrap()
        );

        // - No polynomial multiplication and expansion, e.g. keeping (x + y) * 2 as-is
        let factor = Add::new(vec![x.clone(), y.clone()]).unwrap();
        let mul6 =
            Mul::new(vec![c1.clone(), factor.clone(), factor.clone(), factor.clone()]).unwrap();
        mul = downcast_from_arc::<Mul>(&mul6).unwrap();

        assert_eq!(mul.factors(), vec![Power::new(factor.clone(), 3).unwrap()]);
    }

    #[test]
    fn test_serialization() {
        let op = Mul::new(vec![
            make_number_i64(256u32),
            make_number_f64(64u32),
            make_number_complex(64u32),
            make_number_rational(256u32),
            make_symbol(4u32),
            Power::new(make_symbol(4u32), rand::random_range(-256..=256) as i64).unwrap(),
        ])
        .unwrap();
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let c1 = make_number_complex(64u32);
        let s1 = make_symbol(4u32);
        let s2 = make_symbol(4u32);
        let op = Mul::new(vec![c1.clone(), s1.clone(), s2.clone()]).unwrap();

        assert!(is_expr_type::<Mul>(&op));
        assert!(!is_zero_expr(&op));
        assert!(!is_one_expr(&op));

        let c2 = make_number_rational(256u32);
        let op1 = Mul::new(vec![c1.clone(), s1.clone(), s2.clone()]).unwrap();
        let op2 = Mul::new(vec![s2.clone(), s1.clone(), c1.clone()]).unwrap();
        let op3 = Mul::new(vec![c2.clone(), s1.clone(), s2.clone()]).unwrap();
        let op4 = Mul::new(vec![c1.clone(), s1.clone()]).unwrap();

        assert!(Arc::ptr_eq(&op, &op1));
        assert!(Arc::ptr_eq(&op, &op2));
        assert!(!Arc::ptr_eq(&op, &op3));
        assert!(!Arc::ptr_eq(&op, &op4));
    }
}
