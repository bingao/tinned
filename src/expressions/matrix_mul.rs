use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{MatrixAdd, Mul, Number, ZeroOperator};
use crate::internal::{intern_expr, multi_expression_format, multi_expression_hash};
use crate::perturbations::Perturbation;
use crate::public::{
    NumberTolerance, downcast_from_arc, downcast_from_ref, expression_error,
    generic_expression_error, is_expr_type, is_one_expr, is_zero_expr, subtract_exprs,
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
                    if num.is_zero(None) {
                        return Ok(ZeroOperator::new());
                    }
                }
                all_coefficients.push(term.clone());
            } else if is_expr_type::<ZeroOperator>(term) {
                return Ok(ZeroOperator::new());
            } else if let Some(mat_mul) = downcast_from_arc::<MatrixMul>(term) {
                all_coefficients.push(mat_mul.coefficient().clone());
                all_factors.extend_from_slice(mat_mul.factors());
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
            1 if is_one_expr(&coefficient, None) => Ok(all_factors.pop().unwrap()),
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

const DEFAULT_HASH_DELIMITER: &str = ";";
const DEFAULT_FMT_DELIMITER: &str = " * ";

impl_mul_traits!(MatrixMul, false, DEFAULT_HASH_DELIMITER, DEFAULT_FMT_DELIMITER);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::Symbol;
    use crate::expressions::exch_corr_energy::test_utils::make_exch_corr_energy;
    use crate::expressions::number::test_utils::make_number_complex;
    use crate::expressions::one_elec_operator::test_utils::make_one_elec_operator;
    use crate::expressions::symbol::test_utils::make_symbol;
    use crate::expressions::two_elec_operator::test_utils::make_two_elec_operator;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use num_complex::Complex64;

    test_struct_safety!(MatrixMul);

    test_thread_interning!({
        MatrixMul::new(vec![
            Symbol::new("w"),
            Number::from_complex(Complex64::new(0.0, -1.0)),
            make_wfn_parameter("psi"),
            make_two_elec_operator("op(2el)", Some(make_wfn_parameter("phi"))),
            MatrixAdd::new(vec![
                make_wfn_parameter("phi"),
                make_two_elec_operator("op(2el)", Some(make_wfn_parameter("psi"))),
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
        let op_b = make_one_elec_operator("");
        let op_c = make_two_elec_operator("", None);

        let mul1 = MatrixMul::new(vec![
            c1.clone(),
            c2.clone(),
            op_a.clone(),
            MatrixAdd::new(vec![op_a.clone(), op_b.clone()]).unwrap(),
            MatrixAdd::new(vec![op_b.clone(), op_c.clone()]).unwrap(),
            op_c.clone(),
        ])
        .unwrap();

        assert!(is_expr_type::<MatrixMul>(&mul1));

        let mut mul = downcast_from_arc::<MatrixMul>(&mul1).unwrap();
        let expected_coef = Mul::new(vec![c1.clone(), c2.clone()]).unwrap();
        let expected_factors = vec![
            op_a.clone(),
            MatrixAdd::new(vec![op_a.clone(), op_b.clone()]).unwrap(),
            MatrixAdd::new(vec![op_b.clone(), op_c.clone()]).unwrap(),
            op_c.clone(),
        ];

        // - Order of input terms is preserved
        assert_eq!(
            mul,
            &MatrixMul {
                coefficient: expected_coef.clone(),
                factors: expected_factors.clone(),
            }
        );
        assert_eq!(mul.coefficient(), &expected_coef);
        assert_eq!(mul.factors(), &expected_factors);

        assert_eq!(
            mul1.hash_key(),
            format!(
                "MatrixMul({}{}{})",
                expected_coef.hash_key(),
                DEFAULT_HASH_DELIMITER,
                multi_expression_hash(&expected_factors, DEFAULT_HASH_DELIMITER),
            )
        );
        assert!(!mul1.is_scalar());

        if is_one_expr(&expected_coef, None) {
            assert_eq!(
                format!("{}", mul1),
                format!("{}", multi_expression_format(&expected_factors, DEFAULT_FMT_DELIMITER))
            );
        } else {
            assert_eq!(
                format!("{}", mul1),
                format!(
                    "{}{}{}",
                    expected_coef,
                    DEFAULT_FMT_DELIMITER,
                    multi_expression_format(&expected_factors, DEFAULT_FMT_DELIMITER),
                )
            );
        }

        let mul2 = MatrixMul::new(vec![
            c1.clone(),
            op_a.clone(),
            MatrixAdd::new(vec![op_a.clone(), op_b.clone()]).unwrap(),
            MatrixAdd::new(vec![op_b.clone(), op_c.clone()]).unwrap(),
            op_c.clone(),
            c2.clone(),
        ])
        .unwrap();
        let mul3 = MatrixMul::new(vec![
            c1.clone(),
            op_a.clone(),
            MatrixMul::new(vec![
                c2.clone(),
                MatrixAdd::new(vec![op_a.clone(), op_b.clone()]).unwrap(),
            ])
            .unwrap(),
            MatrixAdd::new(vec![op_b.clone(), op_c.clone()]).unwrap(),
            op_c.clone(),
        ])
        .unwrap();
        let mul4 = MatrixMul::new(vec![
            op_a.clone(),
            MatrixAdd::new(vec![op_a.clone(), op_b.clone()]).unwrap(),
            MatrixAdd::new(vec![op_b.clone(), op_c.clone()]).unwrap(),
            op_c.clone(),
        ])
        .unwrap();

        assert_eq!(&mul1, &mul2);
        assert_eq!(&mul1, &mul3);
        assert_ne!(&mul1, &mul4);

        mul = downcast_from_arc::<MatrixMul>(&mul4).unwrap();

        assert!(is_one_expr(mul.coefficient(), None));

        // - Remove empty MatrixMul([]) -> op(0)
        assert_eq!(&MatrixMul::new(vec![]).unwrap(), &ZeroOperator::new());

        // - Remove redundant MatrixMul([A]) -> A
        assert_eq!(&MatrixMul::new(vec![op_a.clone()]).unwrap(), &op_a);

        // - Identities, ensure A * 0 = op(0), A * op(0) = op(0)
        assert_eq!(
            &MatrixMul::new(vec![op_a.clone(), Number::zero()]).unwrap(),
            &ZeroOperator::new()
        );
        assert_eq!(
            &MatrixMul::new(vec![op_a.clone(), ZeroOperator::new()]).unwrap(),
            &ZeroOperator::new()
        );

        // - Flatten nested MatrixMul
        // - Numeric simplifications, e.g. 3 * 2 -> 6, (1/2) * 4 -> 2, 3 * 0 -> 0
        // - Series multiplication, e.g. ((2 * A) * (3 * A)) * A -> 6 * A * A * A
        assert_eq!(
            &MatrixMul::new(vec![
                MatrixMul::new(vec![
                    c1.clone(),
                    op_a.clone(),
                    MatrixMul::new(vec![op_b.clone(), op_c.clone()]).unwrap(),
                ])
                .unwrap(),
                MatrixMul::new(vec![
                    c2.clone(),
                    op_a.clone(),
                    MatrixMul::new(vec![op_c.clone(), op_b.clone()]).unwrap(),
                ])
                .unwrap(),
                op_a.clone(),
                op_b.clone(),
            ])
            .unwrap(),
            &MatrixMul::new(vec![
                Mul::new(vec![c1.clone(), c2.clone()]).unwrap(),
                op_a.clone(),
                op_b.clone(),
                op_c.clone(),
                op_a.clone(),
                op_c.clone(),
                op_b.clone(),
                op_a.clone(),
                op_b.clone(),
            ])
            .unwrap()
        );

        // - No polynomial multiplication and expansion, e.g. keeping (A + B) * 2 as-is
        let factor = MatrixAdd::new(vec![op_a.clone(), op_b.clone()]).unwrap();
        let mul5 = MatrixMul::new(vec![c1.clone(), factor.clone(), factor.clone(), factor.clone()])
            .unwrap();
        mul = downcast_from_arc::<MatrixMul>(&mul5).unwrap();

        assert_eq!(mul.factors(), vec![factor.clone(), factor.clone(), factor.clone()]);

        assert_eq!(
            &MatrixMul::new(vec![c1.clone(), c2.clone()]).unwrap(),
            &Mul::new(vec![c1.clone(), c2.clone()]).unwrap()
        );
    }

    #[test]
    fn test_differentiation() {
        let coef = make_exch_corr_energy("", None, None, None);
        let op_a = make_wfn_parameter("");
        let op_b = make_one_elec_operator("");
        let op_c = make_two_elec_operator("", None);
        let mul =
            MatrixMul::new(vec![coef.clone(), op_a.clone(), op_b.clone(), op_c.clone()]).unwrap();

        let p = make_perturbation_symbol(4u32, 4u32);
        let diff_mul = mul.differentiate(&p).unwrap();
        let diff_coef = coef.differentiate(&p).unwrap();
        let diff_a = op_a.differentiate(&p).unwrap();
        let diff_b = op_b.differentiate(&p).unwrap();
        let diff_c = op_c.differentiate(&p).unwrap();

        assert_eq!(
            &diff_mul,
            &MatrixAdd::new(vec![
                MatrixMul::new(vec![diff_coef.clone(), op_a.clone(), op_b.clone(), op_c.clone()])
                    .unwrap(),
                MatrixMul::new(vec![coef.clone(), diff_a.clone(), op_b.clone(), op_c.clone()])
                    .unwrap(),
                MatrixMul::new(vec![coef.clone(), op_a.clone(), diff_b.clone(), op_c.clone()])
                    .unwrap(),
                MatrixMul::new(vec![coef.clone(), op_a.clone(), op_b.clone(), diff_c.clone()])
                    .unwrap(),
            ])
            .unwrap()
        );
    }

    #[test]
    fn test_serialization() {
        let op = MatrixMul::new(vec![
            make_symbol(4u32),
            make_number_complex(64u32),
            make_wfn_parameter(""),
            make_one_elec_operator(""),
            MatrixAdd::new(vec![make_wfn_parameter(""), make_two_elec_operator("", None)]).unwrap(),
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
        let op_b = make_one_elec_operator("");
        let op_c = make_two_elec_operator("", None);

        let mul = MatrixMul::new(vec![
            c1.clone(),
            c2.clone(),
            op_a.clone(),
            MatrixAdd::new(vec![op_a.clone(), op_b.clone()]).unwrap(),
            MatrixAdd::new(vec![op_b.clone(), op_c.clone()]).unwrap(),
            op_c.clone(),
        ])
        .unwrap();

        assert!(is_expr_type::<MatrixMul>(&mul));
        assert!(!is_zero_expr(&mul, None));
        assert!(!is_one_expr(&mul, None));

        let mul1 = MatrixMul::new(vec![
            c1.clone(),
            c2.clone(),
            op_a.clone(),
            MatrixAdd::new(vec![op_a.clone(), op_b.clone()]).unwrap(),
            MatrixAdd::new(vec![op_b.clone(), op_c.clone()]).unwrap(),
            op_c.clone(),
        ])
        .unwrap();
        let mul2 = MatrixMul::new(vec![
            c1.clone(),
            op_a.clone(),
            c2.clone(),
            MatrixAdd::new(vec![op_b.clone(), op_a.clone()]).unwrap(),
            MatrixAdd::new(vec![op_c.clone(), op_b.clone()]).unwrap(),
            op_c.clone(),
        ])
        .unwrap();
        let mul3 = MatrixMul::new(vec![
            c1.clone(),
            op_a.clone(),
            c2.clone(),
            MatrixAdd::new(vec![op_b.clone(), op_a.clone()]).unwrap(),
            op_c.clone(),
            MatrixAdd::new(vec![op_c.clone(), op_b.clone()]).unwrap(),
        ])
        .unwrap();

        assert!(Arc::ptr_eq(&mul, &mul1));
        assert!(Arc::ptr_eq(&mul, &mul2));
        assert!(!Arc::ptr_eq(&mul, &mul3));
    }
}
