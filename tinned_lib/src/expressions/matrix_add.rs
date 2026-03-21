use std::collections::HashMap;
use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{Add, MatrixMul, Number, ZeroOperator};
use crate::internal::{intern_expr, sort_expressions_grouped_by};
use crate::public::{
    downcast_from_arc, expression_error, is_expr_type, is_one_expr, is_zero_expr, unreachable_error,
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
    // - Sort terms based on type names and hash keys
    pub fn new(terms: Vec<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        let mut merged: HashMap<String, (Arc<dyn Expr>, Vec<Arc<dyn Expr>>)> = HashMap::new();

        #[inline]
        fn collect_terms(
            expr: &Arc<dyn Expr>,
            merged: &mut HashMap<String, (Arc<dyn Expr>, Vec<Arc<dyn Expr>>)>,
        ) -> Result<(), TinnedError> {
            if expr.is_scalar() {
                return Err(expression_error(
                    "MatrixAdd::new::collect_terms() get a scalar expr",
                    expr,
                    None,
                ));
            }

            if is_expr_type::<ZeroOperator>(expr) {
                return Ok(());
            } else if let Some(mat_add) = downcast_from_arc::<MatrixAdd>(expr) {
                for term in mat_add.terms() {
                    collect_terms(term, merged)?;
                }
            } else if let Some(mat_mul) = downcast_from_arc::<MatrixMul>(expr) {
                if mat_mul.factors().is_empty() {
                    return Err(unreachable_error(
                        "MatrixAdd::new() got MatrixMul with empty factors",
                        expr,
                        None,
                    ));
                }
                let base_expr = if mat_mul.factors().len() == 1 {
                    mat_mul.factors()[0].clone()
                } else {
                    MatrixMul::new(mat_mul.factors().to_vec())?
                };
                let key = base_expr.hash_key();
                let coef = mat_mul.coefficient();
                if let Some((existing_expr, existing_coef)) = merged.get_mut(&key) {
                    if existing_expr == &base_expr {
                        existing_coef.push(coef.clone());
                        return Ok(());
                    }
                }
                merged.insert(key, (base_expr, vec![coef.clone()]));
            } else {
                let key = expr.hash_key();
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

        let mut simplified_terms: Vec<Arc<dyn Expr>> = Vec::with_capacity(merged.len());

        for (expr, all_coef) in merged.into_values() {
            let coef = Add::new(all_coef)?;
            if !is_zero_expr(&coef, None) {
                if is_one_expr(&coef, None) {
                    simplified_terms.push(expr);
                } else {
                    simplified_terms.push(MatrixMul::new(vec![coef, expr])?);
                }
            }
        }

        let mut sorted_terms = sort_expressions_grouped_by(&simplified_terms, |e| e.type_name());

        match sorted_terms.len() {
            0 => Ok(ZeroOperator::new()),
            1 => Ok(sorted_terms.pop().unwrap()),
            _ => Ok(intern_expr(Arc::new(Self {
                terms: sorted_terms,
            }))),
        }
    }

    #[inline]
    pub fn terms(&self) -> &[Arc<dyn Expr>] {
        &self.terms
    }
}

const DEFAULT_HASH_DELIMITER: &str = ";";
const DEFAULT_FMT_DELIMITER: &str = " + ";

impl_add_traits!(MatrixAdd, DEFAULT_HASH_DELIMITER, DEFAULT_FMT_DELIMITER, false);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::Symbol;
    use crate::expressions::ao_two_elec_matrix::test_utils::make_ao_two_elec_matrix;
    use crate::expressions::number::test_utils::make_number_complex;
    use crate::expressions::one_elec_matrix::test_utils::make_one_elec_matrix;
    use crate::expressions::symbol::test_utils::make_symbol;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::internal::join_mapped;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use num_complex::Complex64;

    test_struct_safety!(MatrixAdd);

    test_thread_interning!({
        MatrixAdd::new(vec![
            MatrixMul::new(vec![
                Number::from_complex(Complex64::new(0.0, -1.0)),
                make_wfn_parameter("psi"),
            ])
            .unwrap(),
            make_wfn_parameter("phi"),
            MatrixMul::new(vec![
                Symbol::new("w"),
                make_ao_two_elec_matrix("op(2el)", Some(make_wfn_parameter("psi"))),
            ])
            .unwrap(),
        ])
        .unwrap()
    });

    #[test]
    fn test_impl_expr() {
        let c1 = make_number_complex(64u32);
        let c2 = make_symbol(4u32);
        let op_a = make_wfn_parameter("");
        let op_b = make_one_elec_matrix("", false);
        let op_c = make_ao_two_elec_matrix("", None);

        let add1 = MatrixAdd::new(vec![
            MatrixMul::new(vec![c1.clone(), op_a.clone()]).unwrap(),
            MatrixMul::new(vec![c2.clone(), op_b.clone()]).unwrap(),
            op_c.clone(),
        ])
        .unwrap();

        assert!(is_expr_type::<MatrixAdd>(&add1));

        let add = downcast_from_arc::<MatrixAdd>(&add1).unwrap();
        let expected_terms = sort_expressions_grouped_by(
            &vec![
                MatrixMul::new(vec![c1.clone(), op_a.clone()]).unwrap(),
                MatrixMul::new(vec![c2.clone(), op_b.clone()]).unwrap(),
                op_c.clone(),
            ],
            |e| e.type_name(),
        );

        // - Sort terms based on type names and hash keys
        assert_eq!(
            add,
            &MatrixAdd {
                terms: expected_terms.clone()
            }
        );
        assert_eq!(add.terms(), &expected_terms);

        assert_eq!(
            add1.hash_key(),
            format!(
                "MatrixAdd({})",
                join_mapped(&expected_terms, DEFAULT_HASH_DELIMITER, |term| term.hash_key())
            )
        );
        assert!(!add1.is_scalar());
        assert_eq!(
            format!("{}", add1),
            format!(
                "({})",
                join_mapped(&expected_terms, DEFAULT_FMT_DELIMITER, |term| term.to_string())
            )
        );

        let add2 = MatrixAdd::new(vec![
            MatrixMul::new(vec![c1.clone(), op_a.clone()]).unwrap(),
            MatrixMul::new(vec![c2.clone(), op_b.clone()]).unwrap(),
            op_c.clone(),
        ])
        .unwrap();
        let add3 = MatrixAdd::new(vec![
            op_a.clone(),
            MatrixMul::new(vec![
                Add::new(vec![c1.clone(), Number::from_i64(-1)]).unwrap(),
                op_a.clone(),
            ])
            .unwrap(),
            op_b.clone(),
            MatrixMul::new(vec![
                Add::new(vec![c2.clone(), Number::from_i64(-1)]).unwrap(),
                op_b.clone(),
            ])
            .unwrap(),
            op_c.clone(),
        ])
        .unwrap();
        let add4 = MatrixAdd::new(vec![
            op_a.clone(),
            MatrixMul::new(vec![c2.clone(), op_b.clone()]).unwrap(),
            op_c.clone(),
        ])
        .unwrap();

        assert_eq!(&add1, &add2);
        assert_eq!(&add1, &add3);
        assert_ne!(&add1, &add4);

        // - Remove empty MatrixAdd([]) -> op(0)
        assert_eq!(&MatrixAdd::new(vec![]).unwrap(), &ZeroOperator::new());

        // - Remove redundant MatrixAdd([A]) -> A
        assert_eq!(&MatrixAdd::new(vec![op_a.clone()]).unwrap(), &op_a);

        // - Identities: A + op(0) = A
        assert_eq!(&MatrixAdd::new(vec![op_a.clone(), ZeroOperator::new()]).unwrap(), &op_a);

        // - Combine like terms: 2*A*B + 3*A*B -> 5*A*B
        let coef3 = make_number_complex(64u32);

        assert_eq!(
            &MatrixAdd::new(vec![
                MatrixMul::new(vec![c1.clone(), op_a.clone(), op_b.clone()]).unwrap(),
                MatrixMul::new(vec![c2.clone(), op_a.clone(), op_b.clone()]).unwrap(),
                MatrixMul::new(vec![coef3.clone(), op_a.clone(), op_b.clone()]).unwrap(),
            ])
            .unwrap(),
            &MatrixMul::new(vec![
                Add::new(vec![c1.clone(), c2.clone(), coef3.clone()]).unwrap(),
                op_a.clone(),
                op_b.clone(),
            ])
            .unwrap()
        );

        // - Flatten nested MatrixAdd: (A + B) + (A + C) -> 2A + B + C
        assert_eq!(
            &MatrixAdd::new(vec![
                MatrixAdd::new(vec![
                    MatrixAdd::new(vec![
                        op_a.clone(),
                        MatrixMul::new(vec![op_b.clone(), op_c.clone()]).unwrap(),
                    ])
                    .unwrap(),
                    MatrixAdd::new(vec![
                        op_a.clone(),
                        MatrixMul::new(vec![op_c.clone(), op_b.clone()]).unwrap(),
                    ])
                    .unwrap(),
                ])
                .unwrap(),
                op_a.clone(),
                op_b.clone(),
            ])
            .unwrap(),
            &MatrixAdd::new(vec![
                MatrixMul::new(vec![Number::from_i64(3), op_a.clone()]).unwrap(),
                MatrixMul::new(vec![op_b.clone(), op_c.clone()]).unwrap(),
                MatrixMul::new(vec![op_c.clone(), op_b.clone()]).unwrap(),
                op_b.clone(),
            ])
            .unwrap()
        );
    }

    #[test]
    fn test_differentiation() {
        let op_a =
            MatrixMul::new(vec![make_number_complex(64u32), make_wfn_parameter("")]).unwrap();
        let op_b =
            MatrixMul::new(vec![make_symbol(4u32), make_one_elec_matrix("", false)]).unwrap();
        let op_c = make_ao_two_elec_matrix("", None);
        let add = MatrixAdd::new(vec![op_a.clone(), op_b.clone(), op_c.clone()]).unwrap();

        let p = make_perturbation_symbol(4u32, 4u32);

        assert_eq!(
            &add.differentiate(&p).unwrap(),
            &MatrixAdd::new(vec![
                op_a.differentiate(&p).unwrap(),
                op_b.differentiate(&p).unwrap(),
                op_c.differentiate(&p).unwrap()
            ])
            .unwrap()
        );
    }

    #[test]
    fn test_serialization() {
        let op = MatrixAdd::new(vec![
            MatrixMul::new(vec![make_number_complex(64u32), make_wfn_parameter("")]).unwrap(),
            make_wfn_parameter(""),
            MatrixMul::new(vec![make_symbol(4u32), make_ao_two_elec_matrix("", None)]).unwrap(),
        ])
        .unwrap();
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let c1 = make_number_complex(64u32);
        let c2 = make_symbol(4u32);
        let op_a = make_wfn_parameter("");
        let op_b = make_one_elec_matrix("", false);
        let op_c = make_ao_two_elec_matrix("", None);

        let add = MatrixAdd::new(vec![
            MatrixMul::new(vec![c1.clone(), op_a.clone()]).unwrap(),
            MatrixMul::new(vec![c2.clone(), op_b.clone()]).unwrap(),
            op_c.clone(),
        ])
        .unwrap();

        assert!(is_expr_type::<MatrixAdd>(&add));
        assert!(!is_zero_expr(&add, None));
        assert!(!is_one_expr(&add, None));

        let add1 = MatrixAdd::new(vec![
            MatrixMul::new(vec![c1.clone(), op_a.clone()]).unwrap(),
            MatrixMul::new(vec![c2.clone(), op_b.clone()]).unwrap(),
            op_c.clone(),
        ])
        .unwrap();
        let add2 = MatrixAdd::new(vec![
            op_c.clone(),
            MatrixMul::new(vec![c2.clone(), op_b.clone()]).unwrap(),
            MatrixMul::new(vec![c1.clone(), op_a.clone()]).unwrap(),
        ])
        .unwrap();
        let add3 = MatrixAdd::new(vec![
            op_a.clone(),
            MatrixMul::new(vec![c2.clone(), op_b.clone()]).unwrap(),
            op_c.clone(),
        ])
        .unwrap();

        assert!(Arc::ptr_eq(&add, &add1));
        assert!(Arc::ptr_eq(&add, &add2));
        assert!(!Arc::ptr_eq(&add, &add3));
    }
}
