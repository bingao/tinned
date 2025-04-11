use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{MatrixMul, OneElecOperator, TemporumOperator, ZeroOperator};
use crate::perturbations::{
    pert_multichain_display, pert_multichain_hash_key, PertMultichain, Perturbation,
};
use crate::utils::{downcast_from_ref, intern_expr, is_expr_type};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct TemporumOverlap {
    braket: Arc<dyn Expr>,
    dependencies: PertMultichain,
    derivative: PertMultichain,
}

impl TemporumOverlap {
    #[inline]
    pub fn builder(dependencies: PertMultichain) -> TemporumOverlapBuilder {
        TemporumOverlapBuilder {
            dependencies,
        }
    }

    #[inline]
    pub fn braket(&self) -> &Arc<dyn Expr> {
        &self.braket
    }

    #[inline]
    pub fn dependencies(&self) -> &PertMultichain {
        &self.dependencies
    }

    #[inline]
    pub fn derivative(&self) -> &PertMultichain {
        &self.derivative
    }
}

#[derive(Debug)]
pub struct TemporumOverlapBuilder {
    dependencies: PertMultichain,
}

impl TemporumOverlapBuilder {
    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        let bra = OneElecOperator::builder("Sb").dependencies(self.dependencies.clone()).build()?;
        let dt_bra = TemporumOperator::builder(bra).on_ket(false).build()?;

        let ket = OneElecOperator::builder("Sk").dependencies(self.dependencies.clone()).build()?;
        let dt_ket = TemporumOperator::builder(ket).on_ket(true).build()?;

        Ok(intern_expr(Arc::new(TemporumOverlap {
            braket: MatrixMul::new(vec![dt_bra, dt_ket])?,
            dependencies: self.dependencies,
            derivative: PertMultichain::new(),
        })))
    }
}

#[typetag::serde]
impl Expr for TemporumOverlap {
    #[inline]
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        // We remove braket here, to be consistent with PartialEq
        format!(
            "TemporumOverlap([{}]; [{}])",
            pert_multichain_hash_key(&self.dependencies),
            pert_multichain_hash_key(&self.derivative),
        )
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        false
    }

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        if let Some(op) = downcast_from_ref::<TemporumOverlap>(other) {
            self == op
        } else {
            false
        }
    }

    #[inline]
    fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }

    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_braket = self.braket.differentiate(s)?;
        if is_expr_type::<ZeroOperator>(&diff_braket) {
            return Ok(diff_braket);
        }

        let mut new_deriv = self.derivative.clone();
        *new_deriv.entry(s.clone()).or_insert(0) += 1;

        Ok(intern_expr(Arc::new(Self {
            braket: diff_braket,
            dependencies: self.dependencies.clone(),
            derivative: new_deriv,
        })))
    }
}

impl PartialEq for TemporumOverlap {
    fn eq(&self, other: &Self) -> bool {
        // We do not need to compare braket
        self.dependencies == other.dependencies && self.derivative == other.derivative
    }
}

impl Eq for TemporumOverlap {}

impl std::fmt::Display for TemporumOverlap {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "op(T)^{}", pert_multichain_display(&self.derivative))
    }
}
