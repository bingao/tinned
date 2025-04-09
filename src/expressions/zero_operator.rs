use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::perturbations::Perturbation;
use crate::utils::intern;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ZeroOperator;

impl ZeroOperator {
    #[inline]
    pub fn new() -> Arc<dyn Expr> {
        intern(Arc::new(Self))
    }
}

impl Expr for ZeroOperator {
    #[inline]
    fn as_any(&self) -> &dyn Any {
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
    fn differentiate(&self, _s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(Self::new())
    }
}

impl Display for ZeroOperator {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "op(0)")
    }
}
