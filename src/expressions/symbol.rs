use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::expressions::Number;
use crate::perturbations::Perturbation;
use crate::utils::intern;

/// A scalar symbolic constant that becomes 0 after differentiation.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Symbol {
    name: String,
}

impl Symbol {
    #[inline]
    pub fn new(name: impl Into<String>) -> Arc<dyn Expr> {
        intern(Arc::new(Self { name: name.into() }))
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }
}

impl Expr for Symbol {
    #[inline]
    fn as_any(&self) -> &dyn Any {
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
    fn differentiate(&self, _s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(0.into())
    }
}

impl Display for Symbol {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "{}", self.name)
    }
}
