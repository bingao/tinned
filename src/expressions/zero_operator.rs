use std::collections::HashSet;
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};

#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ZeroOperator;

impl ZeroOperator {
    #[inline]
    pub fn new() -> Arc<dyn Expr> {
        crate::internal::intern_expr(Arc::new(Self))
    }
}

impl ExprInternal for ZeroOperator {
    #[inline]
    fn clone_expr(&self) -> Arc<dyn Expr> {
        Arc::new(self.clone())
    }

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        other.as_any().downcast_ref::<ZeroOperator>().is_some()
    }

    #[allow(unused_variables)]
    #[inline]
    fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str("op(0)")
    }

    #[inline]
    fn hash_key(&self) -> String {
        "ZeroOperator".to_string()
    }
}

#[typetag::serde]
impl Expr for ZeroOperator {
    impl_expr_common_methods!(false);

    #[allow(unused_variables)]
    #[inline]
    fn differentiate(
        &self,
        _s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(Self::new())
    }

    #[inline]
    fn remove(&self, _set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(self.clone_expr())
    }

    #[inline]
    fn retain(
        &self,
        _set: &HashSet<Arc<dyn Expr>>,
        _exact_equality: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(self.clone_expr())
    }
}

impl std::fmt::Display for ZeroOperator {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str("op(0)")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_struct_safety!(ZeroOperator);

    test_thread_interning!(ZeroOperator::new());

    #[test]
    fn test_impl_expr() {
        let z1 = ZeroOperator::new();

        // Test both `as_any` and `downcast_from_arc`
        let z = downcast_from_arc::<ZeroOperator>(&z1).unwrap();
        assert_eq!(z, &ZeroOperator);

        assert_eq!(z1.hash_key(), "ZeroOperator");
        assert!(!z1.is_scalar());

        assert_eq!(format!("{}", z1), "op(0)");

        let z2 = ZeroOperator::new();
        assert_eq!(&z1, &z2);
    }

    #[test]
    fn test_differentiation() {
        let z = ZeroOperator::new();
        let p = make_perturbation_symbol(4u32, 4u32);
        assert_eq!(&z.differentiate(&p).unwrap(), &ZeroOperator::new());
    }

    #[test]
    fn test_serialization() {
        let z = ZeroOperator::new();
        let json = serde_json::to_string(&z).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&z, &deserialized);
    }

    #[test]
    fn test_utils() {
        let z1 = ZeroOperator::new();

        assert!(is_expr_type::<ZeroOperator>(&z1));
        assert!(is_zero_expr(&z1, None));
        assert!(!is_one_expr(&z1, None));

        let z2 = ZeroOperator::new();

        assert!(Arc::ptr_eq(&z1, &z2));
    }
}
