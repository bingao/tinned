use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::Number;
use crate::utils::{downcast_from_arc, downcast_from_ref, intern_expr, invalid_expression_error};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
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
            0 => Ok(Number::one()),
            1 => Ok(base),
            _ => {
                // Flatten nested powers: (x^a)^b -> x^(a * b)
                if let Some(inner) = downcast_from_arc::<Power>(&base) {
                    let combined_exp = inner.exponent * exponent;
                    return Ok(intern_expr(Arc::new(Self {
                        base: inner.base.clone(),
                        exponent: combined_exp,
                    })));
                }

                Ok(intern_expr(Arc::new(Self { base, exponent })))
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

#[typetag::serde]
impl Expr for Power {
    #[inline]
    fn as_any(&self) -> &dyn std::any::Any {
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

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        if let Some(pow) = downcast_from_ref::<Power>(other) {
            self == pow
        } else {
            false
        }
    }

    #[inline]
    fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }

    fn differentiate(
        &self,
        s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let new_exp = self.exponent - 1;
        let diff_base = self.base.differentiate(s)?;

        crate::expressions::Mul::new(vec![
            Number::from_i64(self.exponent),
            Self::new(self.base.clone(), new_exp)?,
            diff_base,
        ])
    }
}

impl PartialEq for Power {
    fn eq(&self, other: &Self) -> bool {
        self.exponent == other.exponent && &self.base == &other.base
    }
}

impl Eq for Power {}

impl std::fmt::Display for Power {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({})^{}", self.base, self.exponent)
    }
}
