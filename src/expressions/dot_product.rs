use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{HermitianTranspose, Number, Transpose, ZeroOperator};
use crate::utils::is_expr_type;

/// Dot product of a bra and a ket (inner product)
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct DotProduct {
    bra: Arc<dyn Expr>,
    ket: Arc<dyn Expr>,
    is_complex: bool,
}

impl DotProduct {
    pub fn new(
        bra: Arc<dyn Expr>,
        ket: Arc<dyn Expr>,
        is_complex: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if bra.is_scalar() || ket.is_scalar() {
            return Err(crate::utils::invalid_expression_error(
                "DotProduct::new() - both arguments must be non-scalar",
                if bra.is_scalar() {
                    &bra
                } else {
                    &ket
                },
            ));
        }

        if is_expr_type::<ZeroOperator>(&bra) || is_expr_type::<ZeroOperator>(&ket) {
            return Ok(Number::zero());
        }

        let bra = if is_complex {
            HermitianTranspose::new(bra)?
        } else {
            Transpose::new(bra)?
        };

        Ok(crate::utils::intern_expr(Arc::new(Self {
            bra,
            ket,
            is_complex,
        })))
    }

    #[inline]
    pub fn bra(&self) -> &Arc<dyn Expr> {
        &self.bra
    }

    #[inline]
    pub fn ket(&self) -> &Arc<dyn Expr> {
        &self.ket
    }

    #[inline]
    pub fn is_complex(&self) -> bool {
        self.is_complex
    }
}

#[typetag::serde]
impl Expr for DotProduct {
    #[inline]
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        format!("DotProduct({}; {}; {})", self.bra.hash_key(), self.ket.hash_key(), self.is_complex)
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        true
    }

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        if let Some(dot) = crate::utils::downcast_from_ref::<DotProduct>(other) {
            self == dot
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
        let diff_bra = self.bra.differentiate(s)?;
        let diff_ket = self.ket.differentiate(s)?;

        Self::new(diff_bra, diff_ket, self.is_complex)
    }
}

impl PartialEq for DotProduct {
    fn eq(&self, other: &Self) -> bool {
        self.is_complex == other.is_complex && &self.bra == &other.bra && &self.ket == &other.ket
    }
}

impl Eq for DotProduct {}

impl std::fmt::Display for DotProduct {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "<{}, {}>", self.bra, self.ket)
    }
}
