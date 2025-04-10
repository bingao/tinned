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
