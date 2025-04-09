use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};

#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Composition {
    name: String,
    order: u32,
    inner: Arc<dyn Expr>,
}

impl Composition {
    #[inline]
    pub fn new(name: impl Into<String>, order: u32, inner: Arc<dyn Expr>) -> Arc<dyn Expr> {
        crate::utils::intern(Arc::new(Self { name: name.into(), order, inner }))
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
        format!("Composition({}^{}, {})", self.name, self.order, self.inner.hash_key())
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        true
    }

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        if let Some(comp) = crate::utils::downcast_from_ref::<Composition>(other) {
            self.name == comp.name && self.order == comp.order && self.inner == comp.inner
        } else {
            false
        }
    }

    fn differentiate(
        &self,
        s: &crate::perturbations::Perturbation,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_outer = Self::new(self.name.clone(), self.order + 1, self.inner.clone());
        let diff_inner = self.inner.differentiate(s)?;

        crate::expressions::Mul::new(vec![diff_outer, diff_inner])
    }
}

impl std::fmt::Display for Composition {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.order == 0 {
            write!(f, "{}({})", self.name, self.inner)
        } else {
            write!(f, "{}^({})({})", self.name, self.order, self.inner)
        }
    }
}
