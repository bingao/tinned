use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::expressions::{MatrixAdd, WfnParameter, ZeroOperator};
use crate::perturbations::{
    is_sub_multichain, pert_multichain_display, pert_multichain_hash_key, PertMultichain,
    Perturbation,
};
use crate::utils::{intern, invalid_expression_error, is_zero_expr};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct TwoElecOperator {
    name: String,
    density: Arc<dyn Expr>,
    dependencies: PertMultichain,
    derivative: PertMultichain,
}

impl TwoElecOperator {
    #[inline]
    pub fn builder(name: impl Into<String>, density: Arc<dyn Expr>) -> TwoElecOperatorBuilder {
        TwoElecOperatorBuilder {
            name: name.into(),
            density,
            dependencies: PertMultichain::new(),
            derivative: PertMultichain::new(),
        }
    }

    #[inline]
    fn builder_from_density(&self, density: Arc<dyn Expr>) -> TwoElecOperatorBuilder {
        TwoElecOperatorBuilder {
            name: self.name.clone(),
            density,
            dependencies: self.dependencies.clone(),
            derivative: self.derivative.clone(),
        }
    }

    #[inline]
    fn builder_from_derivative(&self, derivative: PertMultichain) -> TwoElecOperatorBuilder {
        TwoElecOperatorBuilder {
            name: self.name.clone(),
            density: self.density.clone(),
            dependencies: self.dependencies.clone(),
            derivative,
        }
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[inline]
    pub fn density(&self) -> &Arc<dyn Expr> {
        &self.density
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
pub struct TwoElecOperatorBuilder {
    name: String,
    density: Arc<dyn Expr>,
    dependencies: PertMultichain,
    derivative: PertMultichain,
}

impl TwoElecOperatorBuilder {
    #[inline]
    pub fn dependencies(mut self, deps: PertMultichain) -> Self {
        self.dependencies = deps;
        self
    }

    #[inline]
    pub fn derivative(mut self, deriv: PertMultichain) -> Self {
        self.derivative = deriv;
        self
    }

    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.density.is::<WfnParameter>() {
            if is_sub_multichain(&self.derivative, &self.dependencies) {
                Ok(intern(Arc::new(TwoElecOperator {
                    name: self.name,
                    density: self.density,
                    dependencies: self.dependencies,
                    derivative: self.derivative,
                })))
            } else {
                Ok(ZeroOperator::new())
            }
        } else {
            Err(invalid_expression_error(
                "TwoElecOperatorBuilder::build() - density must be WfnParameter",
                &self.density,
            ))
        }
    }
}

impl Expr for TwoElecOperator {
    #[inline]
    fn as_any(&self) -> &dyn Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "TwoElecOperator({}; {}; [{}]; [{}])",
            self.name,
            self.density.hash_key(),
            pert_multichain_hash_key(&self.dependencies),
            pert_multichain_hash_key(&self.derivative),
        )
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        false
    }

    fn differentiate(&self, s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_density = self.density.differentiate(s)?;

        let term1 = self.builder_from_density(diff_density).build()?;

        let mut new_deriv = self.derivative.clone();
        *new_deriv.entry(s.clone()).or_insert(0) += 1;

        let term2 = self.builder_from_derivative(new_deriv).build()?;

        if is_zero_expr(&term2) {
            return Ok(term1);
        }

        MatrixAdd::new(vec![term1, term2])
    }
}

impl Display for TwoElecOperator {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "{}^{}[{}]", self.name, pert_multichain_display(&self.derivative), self.density)
    }
}
