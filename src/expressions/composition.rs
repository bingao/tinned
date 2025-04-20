use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Composition {
    name: String,
    order: u32,
    inner: Arc<dyn Expr>,
}

impl Composition {
    #[inline]
    pub fn new(name: impl Into<String>, order: u32, inner: Arc<dyn Expr>) -> Arc<dyn Expr> {
        crate::utils::intern_expr(Arc::new(Self {
            name: name.into(),
            order,
            inner,
        }))
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

#[typetag::serde]
impl Expr for Composition {
    #[inline]
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        format!("Composition({}^{}; {})", self.name, self.order, self.inner.hash_key())
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        true
    }

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        if let Some(comp) = crate::utils::downcast_from_ref::<Composition>(other) {
            self == comp
        } else {
            false
        }
    }

    #[inline]
    fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }

    fn differentiate(
        &self,
        s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_outer = Self::new(self.name.clone(), self.order + 1, self.inner.clone());
        let diff_inner = self.inner.differentiate(s)?;

        crate::expressions::Mul::new(vec![diff_outer, diff_inner])
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
    use crate::expressions::symbol::test_utils::{make_symbol, random_alphanumeric};
    use crate::expressions::Power;
    use crate::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_struct_safety!(Composition);

    test_thread_interning!({
        Composition::new("composition", 1, Power::new(make_symbol(0u32), 4).unwrap())
    });

    #[test]
    fn test_impl_expr() {
        let name = random_alphanumeric(4u32);
        let order = rand::random_range(2..=16) as u32;
        let inner = Power::new(make_symbol(2u32), rand::random_range(2..=16) as i64).unwrap();
        let op1 = Composition::new(name.clone(), order, inner.clone());

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

        let op2 = Composition::new(name.clone(), order, inner.clone());
        let op3 = Composition::new(random_alphanumeric(2u32), order, inner.clone());
        let op4 = Composition::new(name.clone(), order + 1, inner.clone());
        let op5 = Composition::new(
            name.clone(),
            order,
            Power::new(make_symbol(4u32), rand::random_range(2..=16) as i64).unwrap(),
        );

        assert_eq!(&op1, &op2);
        assert_ne!(&op1, &op3);
        assert_ne!(&op1, &op4);
        assert_ne!(&op1, &op5);
    }

    #[test]
    fn test_serialization() {
        let op = Composition::new(
            random_alphanumeric(4u32),
            rand::random_range(2..=16) as u32,
            Power::new(make_symbol(2u32), rand::random_range(2..=16) as i64).unwrap(),
        );
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let name = random_alphanumeric(4u32);
        let order = rand::random_range(2..=16) as u32;
        let inner = Power::new(make_symbol(2u32), rand::random_range(2..=16) as i64).unwrap();
        let op1 = Composition::new(name.clone(), order, inner.clone());

        assert!(is_expr_type::<Composition>(&op1));
        assert!(!is_zero_expr(&op1));
        assert!(!is_one_expr(&op1));

        let op2 = Composition::new(name.clone(), order, inner.clone());
        let op3 = Composition::new(random_alphanumeric(2u32), order, inner.clone());
        let op4 = Composition::new(name.clone(), order + 1, inner.clone());
        let op5 = Composition::new(
            name.clone(),
            order,
            Power::new(make_symbol(4u32), rand::random_range(2..=16) as i64).unwrap(),
        );

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));
        assert!(!Arc::ptr_eq(&op1, &op4));
        assert!(!Arc::ptr_eq(&op1, &op5));
    }
}
