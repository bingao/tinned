use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};

/// A scalar symbolic constant that becomes 0 after differentiation.
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Symbol {
    name: String,
}

impl Symbol {
    #[inline]
    pub fn new(name: impl Into<String>) -> Arc<dyn Expr> {
        crate::utils::intern_expr(Arc::new(Self {
            name: name.into(),
        }))
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }
}

#[typetag::serde]
impl Expr for Symbol {
    #[inline]
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        format!("Symbol({})", self.name)
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        true
    }

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        if let Some(s) = crate::utils::downcast_from_ref::<Symbol>(other) {
            self.name == s.name
        } else {
            false
        }
    }

    #[inline]
    fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }

    #[inline]
    fn differentiate(
        &self,
        _s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(crate::expressions::Number::zero())
    }
}

impl std::fmt::Display for Symbol {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_struct_safety!(Symbol);

    test_thread_interning!(Symbol::new("alpha"));

    // Implementation for Expr
    #[test]
    fn test_impl_expr() {
        let s1 = Symbol::new("alpha");

        let s = downcast_from_arc::<Symbol>(&s1).unwrap();
        assert_eq!(
            s,
            &Symbol {
                name: "alpha".into()
            }
        );
        assert_eq!(s.name(), "alpha");

        assert_eq!(s1.hash_key(), "Symbol(alpha)");
        assert!(s1.is_scalar());

        let s2 = Symbol::new("alpha");
        let s3 = Symbol::new("beta");
        assert!(s1 == s2);
        assert!(s1 != s3);

        assert_eq!(format!("{}", s), "alpha");
    }

    // Test serialization and deserialization via `serde_json`
    #[test]
    fn test_serialization() {
        let s = Symbol::new("alpha");
        let json = serde_json::to_string(&s).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert!(s == deserialized);
    }

    // Test utils
    #[test]
    fn test_utils() {
        let s1 = Symbol::new("alpha");
        let s2 = Symbol::new("alpha");
        let s3 = Symbol::new("beta");

        assert!(Arc::ptr_eq(&s1, &s2));
        assert!(!Arc::ptr_eq(&s1, &s3));

        assert!(is_expr_type::<Symbol>(&s1));
        assert!(!is_zero_expr(&s1));
        assert!(!is_one_expr(&s1));
    }
}
