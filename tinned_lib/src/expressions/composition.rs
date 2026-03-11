use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{Mul, Number};
use crate::perturbations::Perturbation;
use crate::public::{downcast_from_arc, expression_error, generic_expression_error, is_zero_expr};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Composition {
    name: String,
    order: u32,
    inner: Arc<dyn Expr>,
}

impl Composition {
    #[inline]
    pub fn new(
        name: impl Into<String>,
        order: u32,
        inner: Arc<dyn Expr>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if !inner.is_scalar() {
            return Err(expression_error(
                "Composition requires scalar inner function",
                &inner,
                None,
            ));
        }

        if is_zero_expr(&inner, None) {
            return Ok(Number::zero());
        }

        Ok(crate::internal::intern_expr(Arc::new(Self {
            name: name.into(),
            order,
            inner,
        })))
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[inline]
    pub fn order(&self) -> u32 {
        self.order
    }

    #[inline]
    pub fn inner(&self) -> &Arc<dyn Expr> {
        &self.inner
    }
}

impl ExprInternal for Composition {
    impl_unary_expr_internal_methods!(Composition, inner, false, |this: &Composition, arg| {
        Self::new(this.name.clone(), this.order, arg)
    });

    #[inline]
    fn hash_key(&self) -> String {
        format!("Composition({}^{}; {})", self.name, self.order, self.inner.hash_key())
    }

    #[inline]
    fn total_order(&self) -> u32 {
        self.order
    }

    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(comp) = downcast_from_arc::<Composition>(other) {
            self.name == comp.name
                && self.order >= comp.order
                && self.inner.deep_eq_superchains(&comp.inner)
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for Composition {
    impl_unary_expr_common_methods!(Composition, inner, True, |this: &Composition, arg| Self::new(
        this.name.clone(),
        this.order,
        arg
    ));

    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        // Differentiation using the chain rule in calculus
        let diff_outer = Self::new(self.name.clone(), self.order + 1, self.inner.clone())?;
        let diff_inner = self.inner.differentiate(s).map_err(|e| {
            generic_expression_error(
                "Composition::differentiate() failed for inner function",
                self,
                Some(Box::new(e)),
            )
        })?;

        Mul::new(vec![diff_outer, diff_inner])
    }
}

impl PartialEq for Composition {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name && self.order == other.order && &self.inner == &other.inner
    }
}

impl Eq for Composition {}

impl std::fmt::Display for Composition {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.order == 0 {
            write!(f, "{}({})", self.name, self.inner)
        } else {
            write!(f, "{}^({})({})", self.name, self.order, self.inner)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::Power;
    use crate::expressions::exch_corr_energy::test_utils::make_exch_corr_energy;
    use crate::expressions::symbol::test_utils::{make_symbol, random_alphanumeric};
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::{is_expr_type, is_one_expr};

    test_struct_safety!(Composition);

    test_thread_interning!({
        Composition::new("composition", 1, Power::new(make_symbol(0u32), 4).unwrap()).unwrap()
    });

    #[test]
    fn test_impl_expr() {
        let name = random_alphanumeric(4u32);
        let order = rand::random_range(2..=16) as u32;
        let inner = Power::new(make_symbol(2u32), rand::random_range(2..=16) as i64).unwrap();
        let op1 = Composition::new(name.clone(), order, inner.clone()).unwrap();

        let op = downcast_from_arc::<Composition>(&op1).unwrap();
        assert_eq!(
            op,
            &Composition {
                name: name.clone(),
                order,
                inner: inner.clone(),
            }
        );
        assert_eq!(op.name(), &name);
        assert_eq!(op.order(), order);
        assert_eq!(op.inner(), &inner);

        assert_eq!(
            op1.hash_key(),
            format!("Composition({}^{}; {})", name, order, inner.hash_key())
        );
        assert!(op1.is_scalar());
        if order == 0 {
            assert_eq!(format!("{}", op1), format!("{}({})", name, inner));
        } else {
            assert_eq!(format!("{}", op1), format!("{}^({})({})", name, order, inner));
        }

        let op2 = Composition::new(name.clone(), order, inner.clone()).unwrap();
        let op3 = Composition::new(random_alphanumeric(2u32), order, inner.clone()).unwrap();
        let op4 = Composition::new(name.clone(), order + 1, inner.clone()).unwrap();
        let op5 = Composition::new(
            name.clone(),
            order,
            Power::new(make_symbol(4u32), rand::random_range(2..=16) as i64).unwrap(),
        )
        .unwrap();

        assert_eq!(&op1, &op2);
        assert_ne!(&op1, &op3);
        assert_ne!(&op1, &op4);
        assert_ne!(&op1, &op5);
    }

    #[test]
    fn test_differentiation() {
        let name = random_alphanumeric(4u32);
        let order: u32 = rand::random_range(2..=16);
        let base = make_exch_corr_energy("", None, None, None);
        let exponent: i64 = rand::random_range(2..=16);
        let inner = Power::new(base.clone(), exponent).unwrap();
        let op = Composition::new(name.clone(), order, inner.clone()).unwrap();

        let p = make_perturbation_symbol(4u32, 4u32);
        let diff_op = op.differentiate(&p).unwrap();

        assert_eq!(
            &diff_op,
            &Mul::new(vec![
                Composition::new(name.clone(), order + 1, inner.clone()).unwrap(),
                inner.differentiate(&p).unwrap(),
            ])
            .unwrap()
        );
    }

    #[test]
    fn test_serialization() {
        let op = Composition::new(
            random_alphanumeric(4u32),
            rand::random_range(2..=16) as u32,
            Power::new(make_symbol(2u32), rand::random_range(2..=16) as i64).unwrap(),
        )
        .unwrap();
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let name = random_alphanumeric(4u32);
        let order = rand::random_range(2..=16) as u32;
        let inner = Power::new(make_symbol(2u32), rand::random_range(2..=16) as i64).unwrap();
        let op1 = Composition::new(name.clone(), order, inner.clone()).unwrap();

        assert!(is_expr_type::<Composition>(&op1));
        assert!(!is_zero_expr(&op1, None));
        assert!(!is_one_expr(&op1, None));

        let op2 = Composition::new(name.clone(), order, inner.clone()).unwrap();
        let op3 = Composition::new(random_alphanumeric(2u32), order, inner.clone()).unwrap();
        let op4 = Composition::new(name.clone(), order + 1, inner.clone()).unwrap();
        let op5 = Composition::new(
            name.clone(),
            order,
            Power::new(make_symbol(4u32), rand::random_range(2..=16) as i64).unwrap(),
        )
        .unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));
        assert!(!Arc::ptr_eq(&op1, &op4));
        assert!(!Arc::ptr_eq(&op1, &op5));
    }
}
