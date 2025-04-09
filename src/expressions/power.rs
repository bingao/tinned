use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::expressions::{Mul, Number};
use crate::perturbations::Perturbation;
use crate::utils::{downcast_expr, intern, invalid_expression_error};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Power {
    base: Arc<dyn Expr>,
    exponent: i64,
}

impl Power {
    pub fn new(base: Arc<dyn Expr>, exponent: i64) -> Result<Arc<dyn Expr>, TinnedError> {
        if !base.is_scalar() {
            return Err(invalid_expression_error("Power::new() - base must be scalar", &base));
        }

        match exponent {
            0 => Ok(1.into()),
            1 => Ok(base),
            _ => {
                // Flatten nested powers: (x^a)^b -> x^(a * b)
                if let Some(inner) = downcast_expr::<Power>(&base) {
                    let combined_exp = inner.exponent * exponent;
                    return Ok(intern(Arc::new(Self {
                        base: inner.base.clone(),
                        exponent: combined_exp,
                    })));
                }

                Ok(intern(Arc::new(Self { base, exponent })))
            },
        }
    }

    #[inline]
    pub fn base(&self) -> &Arc<dyn Expr> {
        &self.base
    }

    #[inline]
    pub fn exponent(&self) -> i64 {
        self.exponent
    }
}

impl Expr for Power {
    #[inline]
    fn as_any(&self) -> &dyn Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        format!("Power({}; {})", self.base.hash_key(), self.exponent)
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        true
    }

    fn differentiate(&self, s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError> {
        let new_exp = self.exponent - 1;
        let diff_base = self.base.differentiate(s)?;

        Mul::new(vec![self.exponent.into(), Self::new(self.base.clone(), new_exp)?, diff_base])
    }
}

impl Display for Power {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "({})^{}", self.base, self.exponent)
    }
}
