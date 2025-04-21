use std::collections::HashMap;
use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Mul, Number};
use crate::utils::{
    downcast_from_arc, downcast_from_ref, intern_expr, invalid_expression_error,
    join_exprs_for_display, join_exprs_for_hash, unreachable_error,
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

const DEFAULT_HASH_DELIMITER: &str = ";";
const DEFAULT_FMT_DELIMITER: &str = " + ";

impl_add_traits!(Add, DEFAULT_HASH_DELIMITER, DEFAULT_FMT_DELIMITER, true);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::exch_corr_energy::test_utils::make_exch_corr_energy;
    use crate::expressions::number::test_utils::{
        make_number_complex, make_number_f64, make_number_i64, make_number_rational,
    };
    use crate::expressions::symbol::test_utils::make_symbol;
    use crate::expressions::two_elec_energy::test_utils::make_two_elec_energy;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::expressions::{Power, Symbol, Trace};
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::utils::{is_expr_type, is_one_expr, is_zero_expr};
    use num_complex::Complex64;
    use num_rational::Rational64;

    test_struct_safety!(Add);

    test_thread_interning!({
        Add::new(vec![
            Number::one(),
            Number::from_f64(3.14),
            Number::from_complex(Complex64::new(0.0, -1.0)),
            Number::from_rational(Rational64::new(22, 7)),
            Symbol::new("x"),
            Mul::new(vec![Number::from_f64(3.14), Symbol::new("x")]).unwrap(),
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
        let add1 = Add::new(vec![c1.clone(), x.clone(), y.clone(), z.clone()]).unwrap();

        assert!(is_expr_type::<Add>(&add1));

        let c1_cast = downcast_from_arc::<Number>(&c1).unwrap();
        let add = downcast_from_arc::<Add>(&add1).unwrap();
        let mut asc_terms = vec![c1.clone(), x.clone(), y.clone(), z.clone()];
        let mut desc_terms = vec![c1.clone(), x.clone(), y.clone(), z.clone()];

        asc_terms.sort_by_key(|f| f.fast_hash());
        desc_terms.sort_by_key(|f| std::cmp::Reverse(f.fast_hash()));

        // - Sort terms based on hash values
        assert_eq!(
            add,
            &Add {
                terms: asc_terms.clone()
            }
        );
        assert_eq!(add.terms(), &asc_terms);
        assert_ne!(add.terms(), &desc_terms);

        assert_eq!(
            add1.hash_key(),
            format!("Add({})", join_exprs_for_hash(&asc_terms, DEFAULT_HASH_DELIMITER))
        );
        assert!(add1.is_scalar());
        assert_eq!(
            format!("{}", add1),
            format!("({})", join_exprs_for_display(&asc_terms, DEFAULT_FMT_DELIMITER))
        );

        let add2 = Add::new(vec![c1.clone(), x.clone(), y.clone(), z.clone()]).unwrap();
        let add3 = Add::new(vec![x.clone(), y.clone(), z.clone()]).unwrap();
        let add4 = Add::new(vec![c1.clone(), x.clone(), y.clone()]).unwrap();

        assert_eq!(&add1, &add2);
        assert_ne!(&add1, &add3);
        assert_ne!(&add1, &add4);

        // - Remove empty Add([]) -> 0
        assert_eq!(&Add::new(vec![]).unwrap(), &Number::zero());

        // - Remove redundant Add([x]) -> x
        assert_eq!(&Add::new(vec![x.clone()]).unwrap(), &x);

        // - Identities: x + 0 = x
        assert_eq!(&Add::new(vec![x.clone(), Number::zero()]).unwrap(), &x);

        // - Numeric simplifications: 3 + 5 -> 8
        let c2 = make_number_complex(64u32);
        let c3 = make_number_complex(64u32);
        let add5 = Add::new(vec![c1.clone(), c2.clone(), c3.clone()]).unwrap();

        assert!(is_expr_type::<Number>(&add5));

        let num = downcast_from_arc::<Number>(&add5).unwrap();
        let c2_cast = downcast_from_arc::<Number>(&c2).unwrap();
        let c3_cast = downcast_from_arc::<Number>(&c3).unwrap();

        assert_eq!(num, &c1_cast.add(&c2_cast.add(&c3_cast)));

        // - Combine like terms: 2*x*y + 3*x*y -> 5*x*y
        assert_eq!(
            &Add::new(vec![
                Mul::new(vec![c1.clone(), x.clone(), y.clone()]).unwrap(),
                Mul::new(vec![c2.clone(), x.clone(), y.clone()]).unwrap(),
                Mul::new(vec![c3.clone(), x.clone(), y.clone()]).unwrap(),
            ])
            .unwrap(),
            &Mul::new(vec![
                Add::new(vec![c1.clone(), c2.clone(), c3.clone()]).unwrap(),
                x.clone(),
                y.clone(),
            ])
            .unwrap()
        );

        // - Flatten nested Add: (x + 2) + (x + 3) -> 2x + 5
        assert_eq!(
            &Add::new(vec![
                Add::new(vec![
                    Add::new(vec![c1.clone(), x.clone()]).unwrap(),
                    Add::new(vec![c2.clone(), x.clone()]).unwrap(),
                ])
                .unwrap(),
                c3.clone(),
                x.clone()
            ])
            .unwrap(),
            &Add::new(vec![
                c1.clone(),
                c2.clone(),
                c3.clone(),
                Mul::new(vec![x.clone(), Number::from_i64(3)]).unwrap(),
            ])
            .unwrap()
        );
    }

    #[test]
    fn test_differentiation() {
        let op_a = Trace::new(make_wfn_parameter("")).unwrap();
        let op_b = make_two_elec_energy("", None, None);
        let op_c = make_exch_corr_energy("", None, None, None);
        let add = Add::new(vec![op_a.clone(), op_b.clone(), op_c.clone()]).unwrap();

        let p = make_perturbation_symbol(4u32, 4u32);

        assert_eq!(
            &add.differentiate(&p).unwrap(),
            &Add::new(vec![
                op_a.differentiate(&p).unwrap(),
                op_b.differentiate(&p).unwrap(),
                op_c.differentiate(&p).unwrap()
            ])
            .unwrap()
        );
    }

    #[test]
    fn test_serialization() {
        let op = Add::new(vec![
            make_number_i64(256u32),
            make_number_f64(64u32),
            make_number_complex(64u32),
            make_number_rational(256u32),
            make_symbol(4u32),
            Mul::new(vec![make_number_f64(64u32), make_symbol(4u32)]).unwrap(),
            Power::new(make_symbol(4u32), rand::random_range(-256..=256) as i64).unwrap(),
            Mul::new(vec![make_number_complex(64u32), make_symbol(4u32)]).unwrap(),
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
        let op = Add::new(vec![c1.clone(), s1.clone(), s2.clone()]).unwrap();

        assert!(is_expr_type::<Add>(&op));
        assert!(!is_zero_expr(&op));
        assert!(!is_one_expr(&op));

        let c2 = make_number_rational(256u32);
        let op1 = Add::new(vec![c1.clone(), s1.clone(), s2.clone()]).unwrap();
        let op2 = Add::new(vec![s2.clone(), s1.clone(), c1.clone()]).unwrap();
        let op3 = Add::new(vec![c2.clone(), s1.clone(), s2.clone()]).unwrap();
        let op4 = Add::new(vec![c1.clone(), s1.clone()]).unwrap();

        assert!(Arc::ptr_eq(&op, &op1));
        assert!(Arc::ptr_eq(&op, &op2));
        assert!(!Arc::ptr_eq(&op, &op3));
        assert!(!Arc::ptr_eq(&op, &op4));
    }
}
