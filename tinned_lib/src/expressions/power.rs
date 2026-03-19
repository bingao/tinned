use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{Mul, Number};
use crate::internal::intern_expr;
use crate::perturbations::Perturbation;
use crate::public::{
    NumberTolerance, downcast_from_arc, expression_error, generic_expression_error,
};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Power {
    base: Arc<dyn Expr>,
    exponent: i64,
}

impl Power {
    pub fn new(base: Arc<dyn Expr>, exponent: i64) -> Result<Arc<dyn Expr>, TinnedError> {
        if !base.is_scalar() {
            return Err(expression_error("Power::new() - base must be scalar", &base, None));
        }

        match exponent {
            0 => Ok(Number::one()),
            1 => Ok(base),
            _ => {
                if let Some(num) = downcast_from_arc::<Number>(&base) {
                    return Ok(num.pow_i64(exponent)?.into());
                }

                // Flatten nested powers: (x^a)^b -> x^(a * b)
                if let Some(inner) = downcast_from_arc::<Power>(&base) {
                    let combined_exp = inner.exponent * exponent;
                    return Ok(intern_expr(Arc::new(Self {
                        base: inner.base.clone(),
                        exponent: combined_exp,
                    })));
                }

                // If `base` is a `Mul`, handle its coefficient separately
                if let Some(mul) = downcast_from_arc::<Mul>(&base) {
                    if !mul.coefficient().is_one(None) {
                        let coeff_power = mul.coefficient().pow_i64(exponent)?.into();
                        let factors_power = intern_expr(Arc::new(Self {
                            base: Mul::new(mul.factors().to_vec())?,
                            exponent,
                        }));
                        return Mul::new(vec![coeff_power, factors_power]);
                    }
                }

                Ok(intern_expr(Arc::new(Self {
                    base,
                    exponent,
                })))
            },
        }
    }

    #[inline]
    pub fn base(&self) -> &Arc<dyn Expr> {
        &self.base
    }

    #[inline]
    pub fn exponent(&self) -> i64 {
        self.exponent
    }
}

impl ExprInternal for Power {
    impl_unary_expr_internal_methods!(Power, base, false, |this: &Power, arg| Self::new(
        arg,
        this.exponent
    ));

    #[inline]
    fn hash_key(&self) -> String {
        format!("Power({}; {})", self.base.hash_key(), self.exponent)
    }

    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(pow) = downcast_from_arc::<Power>(other) {
            self.exponent == pow.exponent && self.base.deep_eq_superchains(&pow.base)
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for Power {
    impl_unary_expr_common_methods!(Power, base, True, |this: &Power, arg| Self::new(
        arg,
        this.exponent
    ));

    #[inline]
    fn apply_zero_rules(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        impl_unary_expr_arg_operation!(
            self,
            base,
            |arg: &Arc<dyn Expr>| arg.apply_zero_rules(freq_tol),
            "Power::apply_zero_rules() failed",
            |this: &Power, arg| Self::new(arg, this.exponent)
        )
    }

    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        let new_exp = self.exponent - 1;
        let diff_base = self.base.differentiate(s).map_err(|e| {
            generic_expression_error(
                "Power::differentiate() failed for base",
                self,
                Some(Box::new(e)),
            )
        })?;

        crate::expressions::Mul::new(vec![
            Number::from_i64(self.exponent),
            Self::new(self.base.clone(), new_exp)?,
            diff_base,
        ])
    }
}

impl PartialEq for Power {
    fn eq(&self, other: &Self) -> bool {
        self.exponent == other.exponent && &self.base == &other.base
    }
}

impl Eq for Power {}

impl std::fmt::Display for Power {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({})^{}", self.base, self.exponent)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::exch_corr_energy::test_utils::make_exch_corr_energy;
    use crate::expressions::symbol::test_utils::make_symbol;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::{is_expr_type, is_one_expr, is_zero_expr};

    test_struct_safety!(Power);

    test_thread_interning!({ Power::new(make_symbol(0u32), 4).unwrap() });

    #[test]
    fn test_impl_expr() {
        let x1 = make_symbol(2u32);
        let op1 = Power::new(x1.clone(), 0).unwrap();
        assert!(is_one_expr(&op1, None));

        let op2 = Power::new(x1.clone(), 1).unwrap();
        assert_eq!(&op2, &x1);

        let exponent1: i64 = rand::random_range(2..=16);
        let op3 = Power::new(x1.clone(), exponent1).unwrap();

        let op = downcast_from_arc::<Power>(&op3).unwrap();
        assert_eq!(
            op,
            &Power {
                base: x1.clone(),
                exponent: exponent1,
            }
        );
        assert_eq!(op.base(), &x1);
        assert_eq!(op.exponent(), exponent1);

        assert_eq!(op3.hash_key(), format!("Power({}; {})", x1.hash_key(), exponent1));
        assert!(op3.is_scalar());
        assert_eq!(format!("{}", op3), format!("({})^{}", x1, exponent1));

        let x2 = make_symbol(4u32);
        let exponent2: i64 = rand::random_range(-32..=32);

        let op4 = Power::new(x1.clone(), exponent1).unwrap();
        let op5 = Power::new(x2.clone(), exponent1).unwrap();
        let op6 = Power::new(x2.clone(), exponent2).unwrap();

        assert_eq!(&op3, &op4);
        assert_ne!(&op3, &op5);
        assert_ne!(&op3, &op6);

        let op7 = Power::new(op3.clone(), exponent2).unwrap();

        assert_eq!(&op7, &Power::new(x1.clone(), exponent1 * exponent2).unwrap());
    }

    #[test]
    fn test_differentiation() {
        let mut op = Power::new(make_symbol(2u32), rand::random_range(2..=16) as i64).unwrap();
        let p = make_perturbation_symbol(4u32, 4u32);

        assert!(is_zero_expr(&op.differentiate(&p).unwrap(), None));

        let base = make_exch_corr_energy("", None, None, None);
        let mut exponent: i64 = rand::random_range(2..=16);
        op = Power::new(base.clone(), exponent).unwrap();

        assert_eq!(
            &op.differentiate(&p).unwrap(),
            &crate::expressions::Mul::new(vec![
                Number::from_i64(exponent),
                Power::new(base.clone(), exponent - 1).unwrap(),
                base.differentiate(&p).unwrap(),
            ])
            .unwrap()
        );

        exponent = rand::random_range(-16..=-1);
        op = Power::new(base.clone(), exponent).unwrap();

        assert_eq!(
            &op.differentiate(&p).unwrap(),
            &crate::expressions::Mul::new(vec![
                Number::from_i64(exponent),
                Power::new(base.clone(), exponent - 1).unwrap(),
                base.differentiate(&p).unwrap(),
            ])
            .unwrap()
        );
    }

    #[test]
    fn test_serialization() {
        let op = Power::new(make_symbol(2u32), rand::random_range(2..=16) as i64).unwrap();
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let x1 = make_symbol(2u32);
        let x2 = make_symbol(4u32);
        let exponent: i64 = rand::random_range(2..=16);
        let op1 = Power::new(x1.clone(), exponent).unwrap();

        assert!(is_expr_type::<Power>(&op1));
        assert!(!is_zero_expr(&op1, None));
        assert!(!is_one_expr(&op1, None));

        let op2 = Power::new(x1.clone(), exponent).unwrap();
        let op3 = Power::new(x2.clone(), exponent).unwrap();
        let op4 = Power::new(x1.clone(), -exponent).unwrap();
        let op5 = Power::new(x2.clone(), -exponent).unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));
        assert!(!Arc::ptr_eq(&op1, &op4));
        assert!(!Arc::ptr_eq(&op1, &op5));
    }
}
