use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};

#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ZeroOperator;

impl ZeroOperator {
    #[inline]
    pub fn new() -> Arc<dyn Expr> {
        crate::utils::intern_expr(Arc::new(Self))
    }
}

#[typetag::serde]
impl Expr for ZeroOperator {
    #[inline]
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        "ZeroOperator".to_string()
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        false
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

    #[allow(unused_variables)]
    #[inline]
    fn differentiate(
        &self,
        _s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(Self::new())
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
    use crate::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_struct_safety!(ZeroOperator);

    test_thread_interning!(ZeroOperator::new());

    #[test]
    fn test_struct() {
        let z = ZeroOperator;

        assert_eq!(z.hash_key(), "ZeroOperator");
        assert!(!z.is_scalar());
        assert_eq!(format!("{}", z), "op(0)");

        assert_eq!(z, ZeroOperator);
    }

    #[test]
    fn test_impl_expr() {
        let z = ZeroOperator::new();

        assert_eq!(z.hash_key(), "ZeroOperator");
        assert!(!z.is_scalar());

        assert_eq!(format!("{}", z), "op(0)");

        let z1 = ZeroOperator::new();
        assert!(z == z1);
    }

    #[test]
    fn test_serialization() {
        let z = ZeroOperator::new();
        let json = serde_json::to_string(&z).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert!(z == deserialized);
    }

    #[test]
    fn test_utils() {
        let z1 = ZeroOperator::new();
        let z2 = ZeroOperator::new();

        assert!(Arc::ptr_eq(&z1, &z2)); // Interning check

        let z = downcast_from_arc::<ZeroOperator>(&z1).unwrap();
        assert_eq!(z, &ZeroOperator);

        assert!(is_expr_type::<ZeroOperator>(&z1));

        assert!(is_zero_expr(&z1));
        assert!(!is_one_expr(&z1));
    }
}
