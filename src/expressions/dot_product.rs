use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{HermitianTranspose, Transpose, ZeroOperator};

/// Dot product of a bra and a ket (inner product)
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
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
                if bra.is_scalar() { &bra } else { &ket },
            ));
        }

        if bra.is::<ZeroOperator>() || ket.is::<ZeroOperator>() {
            return Ok(0.into());
        }

        let bra = if is_complex { HermitianTranspose::new(bra)? } else { Transpose::new(bra)? };

        Ok(crate::utils::intern(Arc::new(Self { bra, ket, is_complex })))
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
            self.is_complex == dot.is_complex && self.bra == dot.bra && self.ket == dot.ket
        } else {
            false
        }
    }

    fn differentiate(
        &self,
        s: &crate::perturbations::Perturbation,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_bra = self.bra.differentiate(s)?;
        let diff_ket = self.ket.differentiate(s)?;

        Self::new(diff_bra, diff_ket, self.is_complex)
    }
}

impl std::fmt::Display for DotProduct {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "<{}, {}>", self.bra, self.ket)
    }
}
