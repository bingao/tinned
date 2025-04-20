use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Add, Number, WfnParameter};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::utils::{
    downcast_from_ref, intern_expr, invalid_expression_error, is_expr_type, is_zero_expr,
};

/// allow_density_swap means we allow inner_density and outer_density to be
/// interchanged when comparing two TwoElecEnergy instances
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
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
            outer_density: Some(self.outer_density.clone()),
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
            outer_density: Some(outer_density),
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
            outer_density: Some(self.outer_density.clone()),
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

        if !is_expr_type::<WfnParameter>(&self.inner_density) {
            return Err(invalid_expression_error(
                "TwoElecEnergyBuilder::build() - inner_density must be WfnParameter",
                &self.inner_density,
            ));
        }

        if !is_expr_type::<WfnParameter>(&outer) {
            return Err(invalid_expression_error(
                "TwoElecEnergyBuilder::build() - outer_density must be WfnParameter",
                &outer,
            ));
        }

        if !self.dependencies.is_subchain(&self.derivative) {
            return Ok(Number::zero());
        }

        Ok(intern_expr(Arc::new(TwoElecEnergy {
            name: self.name,
            inner_density: self.inner_density,
            outer_density: outer,
            allow_density_swap: self.allow_density_swap,
            dependencies: self.dependencies,
            derivative: self.derivative,
        })))
    }
}

#[typetag::serde]
impl Expr for TwoElecEnergy {
    #[inline]
    fn as_any(&self) -> &dyn std::any::Any {
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
            self.dependencies.hash_key(),
            self.derivative.hash_key(),
        )
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        true
    }

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        if let Some(op) = downcast_from_ref::<TwoElecEnergy>(other) {
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
        new_deriv.insert(s);

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
            &self.inner_density == &other.inner_density
                && &self.outer_density == &other.outer_density
        } else {
            // Accept either order
            (&self.inner_density == &other.inner_density
                && &self.outer_density == &other.outer_density)
                || (&self.inner_density == &other.outer_density
                    && &self.outer_density == &other.inner_density)
        }
    }
}

impl Eq for TwoElecEnergy {}

impl std::fmt::Display for TwoElecEnergy {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
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

        write!(f, "{}^{}[{}; {}]", self.name, self.derivative, inner, outer)
    }
}

#[cfg(test)]
const DEFAULT_OPER_NAME: &str = "E(2el)";

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::pert_multichain::test_utils::{
        make_pert_multichain, make_super_multichain,
    };

    #[inline]
    pub fn make_two_elec_energy(
        name: impl Into<String>,
        inner_density: Option<Arc<dyn Expr>>,
        outer_density: Option<Arc<dyn Expr>>,
    ) -> Arc<dyn Expr> {
        let name: String = name.into();
        let inner = inner_density.unwrap_or_else(|| make_wfn_parameter(""));
        let outer = outer_density.unwrap_or_else(|| inner.clone());
        if name.is_empty() {
            let deriv = make_pert_multichain(2u32, 10u32, 1u32, 10u32);
            let deps = make_super_multichain(&deriv, 1u32);
            TwoElecEnergy::builder(random_alphanumeric(DEFAULT_OPER_NAME.len() as u32 + 1), inner)
                .outer_density(outer)
                .dependencies(deps)
                .derivative(deriv)
                .build()
                .unwrap()
        } else {
            let deriv = make_pert_multichain(0u32, 0u32, 1u32, 0u32);
            let deps = make_super_multichain(&deriv, 1u32);
            TwoElecEnergy::builder(name, inner)
                .outer_density(outer)
                .dependencies(deps)
                .derivative(deriv)
                .build()
                .unwrap()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::pert_multichain::test_utils::{
        make_pert_multichain, make_super_multichain,
    };
    use crate::utils::{downcast_from_arc, is_one_expr};

    test_struct_safety!(TwoElecEnergy);

    test_thread_interning!({
        make_two_elec_energy(DEFAULT_OPER_NAME, Some(make_wfn_parameter("density")), None)
    });

    #[test]
    fn test_impl_expr() {
        let mut inner_density = make_wfn_parameter("");
        let mut outer_density = make_wfn_parameter("");
        if inner_density.hash_key() > outer_density.hash_key() {
            std::mem::swap(&mut inner_density, &mut outer_density);
        }
        let deriv = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
        let op0 = TwoElecEnergy::builder(DEFAULT_OPER_NAME, inner_density.clone())
            .derivative(deriv.clone())
            .build()
            .unwrap();

        assert!(is_zero_expr(&op0));

        let deps = make_super_multichain(&deriv, 1u32);
        let op1 = TwoElecEnergy::builder(DEFAULT_OPER_NAME, inner_density.clone())
            .outer_density(outer_density.clone())
            .allow_density_swap(true)
            .dependencies(deps.clone())
            .derivative(deriv.clone())
            .build()
            .unwrap();

        let op = downcast_from_arc::<TwoElecEnergy>(&op1).unwrap();
        assert_eq!(
            op,
            &TwoElecEnergy {
                name: DEFAULT_OPER_NAME.into(),
                inner_density: inner_density.clone(),
                outer_density: outer_density.clone(),
                allow_density_swap: true,
                dependencies: deps.clone(),
                derivative: deriv.clone()
            }
        );

        assert_eq!(op.name(), DEFAULT_OPER_NAME);
        assert_eq!(op.inner_density(), &inner_density);
        assert_eq!(op.outer_density(), &outer_density);
        assert!(op.allow_density_swap());
        assert_eq!(op.dependencies(), &deps);
        assert_eq!(op.derivative(), &deriv);

        let mut op2 = op.builder_from_inner_density(inner_density.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        op2 = op.builder_from_outer_density(outer_density.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        op2 = op.builder_from_derivative(deriv.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        assert_eq!(
            op1.hash_key(),
            format!(
                "TwoElecEnergy({}; {}; {}; [{}]; [{}])",
                DEFAULT_OPER_NAME,
                inner_density.hash_key(),
                outer_density.hash_key(),
                deps.hash_key(),
                deriv.hash_key()
            )
        );
        assert!(op1.is_scalar());
        assert_eq!(
            format!("{}", op1),
            format!("{}^{}[{}; {}]", DEFAULT_OPER_NAME, deriv, inner_density, outer_density)
        );

        let op3 = TwoElecEnergy::builder(DEFAULT_OPER_NAME, inner_density.clone())
            .outer_density(outer_density.clone())
            .allow_density_swap(true)
            .dependencies(deps.clone())
            .derivative(deriv.clone())
            .build()
            .unwrap();
        let op4 = TwoElecEnergy::builder(
            random_alphanumeric(DEFAULT_OPER_NAME.len() as u32 + 1),
            inner_density.clone(),
        )
        .outer_density(outer_density.clone())
        .allow_density_swap(true)
        .dependencies(deps.clone())
        .derivative(deriv.clone())
        .build()
        .unwrap();
        let op5 = TwoElecEnergy::builder(DEFAULT_OPER_NAME, inner_density.clone())
            .outer_density(outer_density.clone())
            .allow_density_swap(true)
            .dependencies(deps.clone())
            .build()
            .unwrap();
        let op6 = TwoElecEnergy::builder(DEFAULT_OPER_NAME, inner_density.clone())
            .outer_density(outer_density.clone())
            .allow_density_swap(true)
            .dependencies(make_super_multichain(&deriv, 2u32))
            .derivative(deriv.clone())
            .build()
            .unwrap();
        let op7 = TwoElecEnergy::builder(DEFAULT_OPER_NAME, make_wfn_parameter("density"))
            .allow_density_swap(true)
            .dependencies(deps.clone())
            .derivative(deriv.clone())
            .build()
            .unwrap();

        assert_eq!(&op1, &op3);
        assert_ne!(&op1, &op4);
        assert_ne!(&op1, &op5);
        assert_ne!(&op1, &op6);
        assert_ne!(&op1, &op7);
    }

    #[test]
    fn test_serialization() {
        let op =
            make_two_elec_energy("", Some(make_wfn_parameter("")), Some(make_wfn_parameter("")));
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let mut inner_density = make_wfn_parameter("");
        let mut outer_density = make_wfn_parameter("");
        if inner_density.hash_key() > outer_density.hash_key() {
            std::mem::swap(&mut inner_density, &mut outer_density);
        }
        let op1 = make_two_elec_energy(
            DEFAULT_OPER_NAME,
            Some(inner_density.clone()),
            Some(outer_density.clone()),
        );

        assert!(is_expr_type::<TwoElecEnergy>(&op1));
        assert!(!is_zero_expr(&op1));
        assert!(!is_one_expr(&op1));

        let op2 = make_two_elec_energy(
            DEFAULT_OPER_NAME,
            Some(inner_density.clone()),
            Some(outer_density.clone()),
        );
        let op3 = make_two_elec_energy(
            DEFAULT_OPER_NAME,
            Some(outer_density.clone()),
            Some(inner_density.clone()),
        );
        let op4 =
            make_two_elec_energy("", Some(inner_density.clone()), Some(outer_density.clone()));
        let op5 =
            make_two_elec_energy(DEFAULT_OPER_NAME, Some(make_wfn_parameter("density")), None);

        let op = downcast_from_arc::<TwoElecEnergy>(&op1).unwrap();
        let op6 = op
            .builder_from_inner_density(outer_density.clone())
            .outer_density(inner_density.clone())
            .allow_density_swap(false)
            .build()
            .unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(Arc::ptr_eq(&op1, &op3));
        assert!(!Arc::ptr_eq(&op1, &op4));
        assert!(!Arc::ptr_eq(&op1, &op5));
        assert!(!Arc::ptr_eq(&op1, &op6));
        assert_eq!(&op1, &op6);
    }
}
