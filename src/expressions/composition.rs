use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::expressions::Mul;
use crate::perturbations::Perturbation;
use crate::utils::intern;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Composition {
    name: String,
    order: u32,
    inner: Arc<dyn Expr>,
}

impl Composition {
    #[inline]
    pub fn new(name: impl Into<String>, order: u32, inner: Arc<dyn Expr>) -> Arc<dyn Expr> {
        intern(Arc::new(Self { name: name.into(), order, inner }))
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

impl Expr for Composition {
    #[inline]
    fn as_any(&self) -> &dyn Any {
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
        if let Some(comp) = downcast_expr::<Composition>(other) {
            self.name == comp.name && self.order == comp.order && self.inner == comp.inner
        } else {
            false
        }
    }

    fn differentiate(&self, s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_outer = Self::new(self.name.clone(), self.order + 1, self.inner.clone());
        let diff_inner = self.inner.differentiate(s)?;

        Mul::new(vec![diff_outer, diff_inner])
    }
}

impl Display for Composition {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        if self.order == 0 {
            write!(f, "{}({})", self.name, self.inner)
        } else {
            write!(f, "{}^({})({})", self.name, self.order, self.inner)
        }
    }
}
