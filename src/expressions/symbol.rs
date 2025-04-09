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
        crate::utils::intern(Arc::new(Self { name: name.into() }))
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
    fn differentiate(
        &self,
        _s: &crate::perturbations::Perturbation,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(0.into())
    }
}

impl std::fmt::Display for Symbol {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.name)
    }
}
