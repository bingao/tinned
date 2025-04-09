use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::expressions::{Add, Number, WfnParameter};
use crate::perturbations::{
    is_sub_multichain, pert_multichain_display, pert_multichain_hash_key, PertMultichain,
    Perturbation,
};
use crate::utils::{intern, invalid_expression_error, is_zero_expr};

/// allow_density_swap means we allow inner_density and outer_density to be
/// interchanged when comparing two TwoElecEnergy instances
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TwoElecEnergy {
    name: String,
    inner_density: Arc<dyn Expr>,
    outer_density: Arc<dyn Expr>,
    allow_density_swap: bool,
    dependencies: PertMultichain,
    derivative: PertMultichain,
}

impl TwoElecEnergy {
    #[inline]
    pub fn builder(name: impl Into<String>, inner_density: Arc<dyn Expr>) -> TwoElecEnergyBuilder {
        TwoElecEnergyBuilder {
            name: name.into(),
            inner_density,
            outer_density: None,
            allow_density_swap: true,
            dependencies: PertMultichain::new(),
            derivative: PertMultichain::new(),
        }
    }

    #[inline]
    fn builder_from_inner_density(&self, inner_density: Arc<dyn Expr>) -> TwoElecEnergyBuilder {
        TwoElecEnergyBuilder {
            name: self.name.clone(),
            inner_density,
            outer_density: self.outer_density.clone(),
            allow_density_swap: true,
            dependencies: self.dependencies.clone(),
            derivative: self.derivative.clone(),
        }
    }

    #[inline]
    fn builder_from_outer_density(&self, outer_density: Arc<dyn Expr>) -> TwoElecEnergyBuilder {
        TwoElecEnergyBuilder {
            name: self.name.clone(),
            inner_density: self.inner_density.clone(),
            outer_density,
            allow_density_swap: true,
            dependencies: self.dependencies.clone(),
            derivative: self.derivative.clone(),
        }
    }

    #[inline]
    fn builder_from_derivative(&self, derivative: PertMultichain) -> TwoElecEnergyBuilder {
        TwoElecEnergyBuilder {
            name: self.name.clone(),
            inner_density: self.inner_density.clone(),
            outer_density: self.outer_density.clone(),
            allow_density_swap: true,
            dependencies: self.dependencies.clone(),
            derivative,
        }
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[inline]
    pub fn inner_density(&self) -> &Arc<dyn Expr> {
        &self.inner_density
    }

    #[inline]
    pub fn outer_density(&self) -> &Arc<dyn Expr> {
        &self.outer_density
    }

    #[inline]
    pub fn allow_density_swap(&self) -> bool {
        self.allow_density_swap
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
pub struct TwoElecEnergyBuilder {
    name: String,
    inner_density: Arc<dyn Expr>,
    outer_density: Option<Arc<dyn Expr>>,
    allow_density_swap: bool,
    dependencies: PertMultichain,
    derivative: PertMultichain,
}

impl TwoElecEnergyBuilder {
    #[inline]
    pub fn outer_density(mut self, outer: Arc<dyn Expr>) -> Self {
        self.outer_density = Some(outer);
        self
    }

    #[inline]
    pub fn allow_density_swap(mut self, allow_density_swap: bool) -> Self {
        self.allow_density_swap = allow_density_swap;
        self
    }

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
        let outer = self.outer_density.unwrap_or_else(|| self.inner_density.clone());

        if !self.inner_density.is::<WfnParameter>() {
            return Err(invalid_expression_error(
                "TwoElecEnergyBuilder::build() - inner_density must be WfnParameter",
                self.inner_density,
            ));
        }

        if !outer.is::<WfnParameter>() {
            return Err(invalid_expression_error(
                "TwoElecEnergyBuilder::build() - outer_density must be WfnParameter",
                outer,
            ));
        }

        if !is_sub_multichain(&self.derivative, &self.dependencies) {
            return Ok(0.into());
        }

        Ok(intern(Arc::new(TwoElecEnergy {
            name: self.name,
            inner_density: self.inner_density,
            outer_density: outer,
            allow_density_swap: self.allow_density_swap,
            dependencies: self.dependencies,
            derivative: self.derivative,
        })))
    }
}

impl Expr for TwoElecEnergy {
    #[inline]
    fn as_any(&self) -> &dyn Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        let (inner, outer) = if self.allow_density_swap {
            let h1 = self.inner_density.hash_key();
            let h2 = self.outer_density.hash_key();
            if h1 <= h2 {
                (h1, h2)
            } else {
                (h2, h1)
            }
        } else {
            (self.inner_density.hash_key(), self.outer_density.hash_key())
        };

        format!(
            "TwoElecEnergy({}; {}; {}; [{}]; [{}])",
            self.name,
            inner,
            outer,
            pert_multichain_hash_key(&self.dependencies),
            pert_multichain_hash_key(&self.derivative),
        )
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        true
    }

    fn differentiate(&self, s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_inner = self.inner_density.differentiate(s)?;
        let diff_outer = self.outer_density.differentiate(s)?;

        let mut terms = vec![
            self.builder_from_inner_density(diff_inner)
                .allow_density_swap(self.allow_density_swap)
                .build()?,
            self.builder_from_outer_density(diff_outer)
                .allow_density_swap(self.allow_density_swap)
                .build()?,
        ];

        let mut new_deriv = self.derivative.clone();
        *new_deriv.entry(s.clone()).or_insert(0) += 1;

        let diff_oper = self
            .builder_from_derivative(new_deriv)
            .allow_density_swap(self.allow_density_swap)
            .build()?;
        if !is_zero_expr(&diff_oper) {
            terms.push(diff_oper);
        }

        Add::new(terms)
    }
}

impl PartialEq for TwoElecEnergy {
    fn eq(&self, other: &Self) -> bool {
        // Compare all fixed fields
        if self.name != other.name
            || self.dependencies != other.dependencies
            || self.derivative != other.derivative
        {
            return false;
        }

        // Handle density equality based on swap flags
        if !self.allow_density_swap && !other.allow_density_swap {
            // Strict matching only
            self.inner_density == other.inner_density && self.outer_density == other.outer_density
        } else {
            // Accept either order
            (self.inner_density == other.inner_density && self.outer_density == other.outer_density)
                || (self.inner_density == other.outer_density
                    && self.outer_density == other.inner_density)
        }
    }
}

impl Eq for TwoElecEnergy {}

impl Display for TwoElecEnergy {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        let (inner, outer) = if self.allow_density_swap {
            let h1 = self.inner_density.hash_key();
            let h2 = self.outer_density.hash_key();
            if h1 <= h2 {
                (self.inner_density.as_ref(), self.outer_density.as_ref())
            } else {
                (self.outer_density.as_ref(), self.inner_density.as_ref())
            }
        } else {
            (self.inner_density.as_ref(), self.outer_density.as_ref())
        };

        write!(
            f,
            "{}^{}[{}; {}]",
            self.name,
            pert_multichain_display(&self.derivative),
            inner,
            outer,
        )
    }
}
