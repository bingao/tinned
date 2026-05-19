use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{
    Add, AoTwoElecMatrix, Number, ResidueParameter, WfnParameter, ZeroOperator,
};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{
    downcast_from_arc, expression_error, generic_expression_error, is_expr_type, is_zero_expr,
};

/// In the atomic orbital (AO) representation, we have two-electron energy 0.5 *
/// Trace[G^{AO}(D^{AO}) * D^{AO}], see for example equation (78), J. Chem.
/// Phys. 129, 214108 (2008). Here D^{AO} is the AO denisty matrix.
///
/// `inner_density` is used to build G^{AO}(D^{AO}) while `outer_density` is
/// the one contracted with G^{AO}(D^{AO}). `allow_density_swap` means we allow
/// `inner_density` and `outer_density` to be / interchanged when comparing two
/// `AoTwoElecEnergy` instances.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct AoTwoElecEnergy {
    name: String,
    inner_density: Arc<dyn Expr>,
    outer_density: Arc<dyn Expr>,
    allow_density_swap: bool,
    dependencies: PertMultichain,
    derivative: PertMultichain,
}

impl AoTwoElecEnergy {
    #[inline]
    pub fn builder(
        name: impl Into<String>,
        inner_density: Arc<dyn Expr>,
    ) -> AoTwoElecEnergyBuilder {
        AoTwoElecEnergyBuilder {
            name: name.into(),
            inner_density,
            outer_density: None,
            allow_density_swap: true,
            dependencies: PertMultichain::new(),
            derivative: PertMultichain::new(),
        }
    }

    #[inline]
    pub fn builder_from_operator(two_elec_op: &AoTwoElecMatrix) -> AoTwoElecEnergyBuilder {
        AoTwoElecEnergyBuilder {
            name: two_elec_op.name().into(),
            inner_density: two_elec_op.density().clone(),
            outer_density: Some(two_elec_op.density().clone()),
            allow_density_swap: true,
            dependencies: two_elec_op.dependencies().clone(),
            derivative: two_elec_op.derivative().clone(),
        }
    }

    #[inline]
    fn with_inner_density(&self, inner_density: Arc<dyn Expr>) -> AoTwoElecEnergyBuilder {
        AoTwoElecEnergyBuilder {
            name: self.name.clone(),
            inner_density,
            outer_density: Some(self.outer_density.clone()),
            allow_density_swap: true,
            dependencies: self.dependencies.clone(),
            derivative: self.derivative.clone(),
        }
    }

    #[inline]
    fn with_outer_density(&self, outer_density: Arc<dyn Expr>) -> AoTwoElecEnergyBuilder {
        AoTwoElecEnergyBuilder {
            name: self.name.clone(),
            inner_density: self.inner_density.clone(),
            outer_density: Some(outer_density),
            allow_density_swap: true,
            dependencies: self.dependencies.clone(),
            derivative: self.derivative.clone(),
        }
    }

    #[inline]
    fn with_derivative(&self, derivative: PertMultichain) -> AoTwoElecEnergyBuilder {
        AoTwoElecEnergyBuilder {
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

    #[inline]
    fn eq_density(&self, other: &Self) -> bool {
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

#[derive(Debug)]
pub struct AoTwoElecEnergyBuilder {
    name: String,
    inner_density: Arc<dyn Expr>,
    outer_density: Option<Arc<dyn Expr>>,
    allow_density_swap: bool,
    dependencies: PertMultichain,
    derivative: PertMultichain,
}

impl AoTwoElecEnergyBuilder {
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
    pub fn derivative(mut self, derivative: PertMultichain) -> Self {
        self.derivative = derivative;
        self
    }

    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        let outer = self.outer_density.unwrap_or_else(|| self.inner_density.clone());

        if is_expr_type::<ZeroOperator>(&self.inner_density) || is_expr_type::<ZeroOperator>(&outer)
        {
            return Ok(Number::zero());
        }

        if !is_expr_type::<WfnParameter>(&self.inner_density)
            && !is_expr_type::<ResidueParameter>(&self.inner_density)
        {
            return Err(expression_error(
                "AoTwoElecEnergyBuilder::build() - inner density must be WfnParameter or ResidueParameter",
                &self.inner_density,
                None,
            ));
        }

        if !is_expr_type::<WfnParameter>(&outer) && !is_expr_type::<ResidueParameter>(&outer) {
            return Err(expression_error(
                "AoTwoElecEnergyBuilder::build() - outer density must be WfnParameter or ResidueParameter",
                &outer,
                None,
            ));
        }

        if !self.dependencies.is_subchain(&self.derivative) {
            return Ok(Number::zero());
        }

        Ok(crate::internal::intern_expr(Arc::new(AoTwoElecEnergy {
            name: self.name,
            inner_density: self.inner_density,
            outer_density: outer,
            allow_density_swap: self.allow_density_swap,
            dependencies: self.dependencies,
            derivative: self.derivative,
        })))
    }
}

impl ExprInternal for AoTwoElecEnergy {
    impl_binary_expr_internal_methods!(
        AoTwoElecEnergy,
        true,
        inner_density,
        outer_density,
        true,
        |this: &AoTwoElecEnergy, inner_density, outer_density| this
            .with_inner_density(inner_density)
            .outer_density(outer_density)
            .allow_density_swap(this.allow_density_swap)
            .build()
    );

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
            "AoTwoElecEnergy({}; {}; {}; [{}]; [{}])",
            self.name,
            inner,
            outer,
            self.dependencies.hash_key(),
            self.derivative.hash_key(),
        )
    }

    #[inline]
    fn expr_order(&self) -> u32 {
        self.derivative.total_order()
    }

    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<AoTwoElecEnergy>(other) {
            // Compare all fixed fields
            if self.name != op.name
                || self.dependencies != op.dependencies
                || !self.derivative.is_subchain(&op.derivative)
            {
                return false;
            }

            // Handle density equality based on swap flags
            if !self.allow_density_swap && !op.allow_density_swap {
                // Strict matching only
                self.inner_density.deep_eq_superchains(&op.inner_density)
                    && self.outer_density.deep_eq_superchains(&op.outer_density)
            } else {
                // Accept either order
                (self.inner_density.deep_eq_superchains(&op.inner_density)
                    && self.outer_density.deep_eq_superchains(&op.outer_density))
                    || (self.inner_density.deep_eq_superchains(&op.outer_density)
                        && self.outer_density.deep_eq_superchains(&op.inner_density))
            }
        } else {
            false
        }
    }

    #[inline]
    fn eq_by_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        // For unambiguous replacement, we require equality of density
        // matrices, and make replacement by considering only derivative of
        // electron repulsion integrals (ERIs).
        if let Some(op) = downcast_from_arc::<AoTwoElecEnergy>(other) {
            self.name == op.name
                && self.dependencies == op.dependencies
                && self.derivative.is_subchain(&op.derivative)
                && self.eq_density(op)
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for AoTwoElecEnergy {
    impl_binary_expr_common_methods!(
        AoTwoElecEnergy,
        true,
        inner_density,
        outer_density,
        |this: &AoTwoElecEnergy, inner_density, outer_density| this
            .with_inner_density(inner_density)
            .outer_density(outer_density)
            .allow_density_swap(this.allow_density_swap)
            .build(),
        false
    );

    fn differentiate(&self, s: Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_inner = self.inner_density.differentiate(s.clone()).map_err(|e| {
            generic_expression_error(
                "AoTwoElecEnergy::differentiate() failed for inner density",
                self,
                Some(Box::new(e)),
            )
        })?;
        let diff_outer = self.outer_density.differentiate(s.clone()).map_err(|e| {
            generic_expression_error(
                "AoTwoElecEnergy::differentiate() failed for outer density",
                self,
                Some(Box::new(e)),
            )
        })?;

        let mut terms = vec![
            self.with_inner_density(diff_inner)
                .allow_density_swap(self.allow_density_swap)
                .build()?,
            self.with_outer_density(diff_outer)
                .allow_density_swap(self.allow_density_swap)
                .build()?,
        ];

        let new_deriv = self.derivative.with_added_perturbation(s);
        let diff_oper =
            self.with_derivative(new_deriv).allow_density_swap(self.allow_density_swap).build()?;

        if !is_zero_expr(&diff_oper, None) {
            terms.push(diff_oper);
        }

        Add::new(terms)
    }
}

impl PartialEq for AoTwoElecEnergy {
    fn eq(&self, other: &Self) -> bool {
        // Compare all fixed fields
        if self.name != other.name
            || self.dependencies != other.dependencies
            || self.derivative != other.derivative
        {
            return false;
        }

        self.eq_density(other)
    }
}

impl Eq for AoTwoElecEnergy {}

impl std::fmt::Display for AoTwoElecEnergy {
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
pub mod test_utils {
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::pert_multichain::test_utils::{
        make_pert_multichain, make_super_multichain,
    };

    pub const TEST_OPER_NAME: &str = "J+K";

    #[inline]
    pub fn make_ao_two_elec_energy(
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
            AoTwoElecEnergy::builder(random_alphanumeric(TEST_OPER_NAME.len() as u32 + 1), inner)
                .outer_density(outer)
                .dependencies(deps)
                .derivative(deriv)
                .build()
                .unwrap()
        } else {
            let deriv = make_pert_multichain(0u32, 0u32, 1u32, 0u32);
            let deps = make_super_multichain(&deriv, 1u32);
            AoTwoElecEnergy::builder(name, inner)
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
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::is_one_expr;

    test_struct_safety!(AoTwoElecEnergy);

    test_thread_interning!({
        make_ao_two_elec_energy(TEST_OPER_NAME, Some(make_wfn_parameter("density")), None)
    });

    #[test]
    fn test_impl_expr() {
        let mut inner_density = make_wfn_parameter("");
        let mut outer_density = make_wfn_parameter("");
        if inner_density.hash_key() > outer_density.hash_key() {
            std::mem::swap(&mut inner_density, &mut outer_density);
        }
        let deriv = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
        let op0 = AoTwoElecEnergy::builder(TEST_OPER_NAME, inner_density.clone())
            .derivative(deriv.clone())
            .build()
            .unwrap();

        assert!(is_zero_expr(&op0, None));

        let deps = make_super_multichain(&deriv, 1u32);
        let op1 = AoTwoElecEnergy::builder(TEST_OPER_NAME, inner_density.clone())
            .outer_density(outer_density.clone())
            .allow_density_swap(true)
            .dependencies(deps.clone())
            .derivative(deriv.clone())
            .build()
            .unwrap();

        let op = downcast_from_arc::<AoTwoElecEnergy>(&op1).unwrap();
        assert_eq!(
            op,
            &AoTwoElecEnergy {
                name: TEST_OPER_NAME.into(),
                inner_density: inner_density.clone(),
                outer_density: outer_density.clone(),
                allow_density_swap: true,
                dependencies: deps.clone(),
                derivative: deriv.clone()
            }
        );

        assert_eq!(op.name(), TEST_OPER_NAME);
        assert_eq!(op.inner_density(), &inner_density);
        assert_eq!(op.outer_density(), &outer_density);
        assert!(op.allow_density_swap());
        assert_eq!(op.dependencies(), &deps);
        assert_eq!(op.derivative(), &deriv);

        let mut op2 = op.with_inner_density(inner_density.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        op2 = op.with_outer_density(outer_density.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        op2 = op.with_derivative(deriv.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        assert_eq!(
            op1.hash_key(),
            format!(
                "AoTwoElecEnergy({}; {}; {}; [{}]; [{}])",
                TEST_OPER_NAME,
                inner_density.hash_key(),
                outer_density.hash_key(),
                deps.hash_key(),
                deriv.hash_key(),
            )
        );
        assert!(op1.is_scalar());
        assert_eq!(
            format!("{}", op1),
            format!("{}^{}[{}; {}]", TEST_OPER_NAME, deriv, inner_density, outer_density)
        );

        let op3 = AoTwoElecEnergy::builder(TEST_OPER_NAME, inner_density.clone())
            .outer_density(outer_density.clone())
            .allow_density_swap(true)
            .dependencies(deps.clone())
            .derivative(deriv.clone())
            .build()
            .unwrap();
        let op4 = AoTwoElecEnergy::builder(
            random_alphanumeric(TEST_OPER_NAME.len() as u32 + 1),
            inner_density.clone(),
        )
        .outer_density(outer_density.clone())
        .allow_density_swap(true)
        .dependencies(deps.clone())
        .derivative(deriv.clone())
        .build()
        .unwrap();
        let op5 = AoTwoElecEnergy::builder(TEST_OPER_NAME, inner_density.clone())
            .outer_density(outer_density.clone())
            .allow_density_swap(true)
            .dependencies(deps.clone())
            .build()
            .unwrap();
        let op6 = AoTwoElecEnergy::builder(TEST_OPER_NAME, inner_density.clone())
            .outer_density(outer_density.clone())
            .allow_density_swap(true)
            .dependencies(make_super_multichain(&deriv, 2u32))
            .derivative(deriv.clone())
            .build()
            .unwrap();
        let op7 = AoTwoElecEnergy::builder(TEST_OPER_NAME, make_wfn_parameter("density"))
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
    fn test_differentiation() {
        let inner_density = make_wfn_parameter("");
        let outer_density = make_wfn_parameter("");
        let len_pert_name: u32 = 2;
        let deps = make_pert_multichain(len_pert_name, 8u32, 1u32, 10u32);
        let op = AoTwoElecEnergy::builder(TEST_OPER_NAME, inner_density.clone())
            .outer_density(outer_density.clone())
            .dependencies(deps.clone())
            .build()
            .unwrap();

        let mut p: Arc<Perturbation> = deps.keys().first().cloned().unwrap();
        let mut diff_op = op.differentiate(p.clone()).unwrap();
        let mut deriv = PertMultichain::new();
        deriv.insert(p.clone());

        assert_eq!(
            &diff_op,
            &Add::new(vec![
                AoTwoElecEnergy::builder(TEST_OPER_NAME, inner_density.clone())
                    .outer_density(outer_density.clone())
                    .dependencies(deps.clone())
                    .derivative(deriv.clone())
                    .build()
                    .unwrap(),
                AoTwoElecEnergy::builder(
                    TEST_OPER_NAME,
                    inner_density.differentiate(p.clone()).unwrap()
                )
                .outer_density(outer_density.clone())
                .dependencies(deps.clone())
                .build()
                .unwrap(),
                AoTwoElecEnergy::builder(TEST_OPER_NAME, inner_density.clone())
                    .outer_density(outer_density.differentiate(p.clone()).unwrap())
                    .dependencies(deps.clone())
                    .build()
                    .unwrap(),
            ])
            .unwrap()
        );

        p = make_perturbation_symbol(len_pert_name + 1u32, 4u32);
        diff_op = op.differentiate(p.clone()).unwrap();

        assert_eq!(
            &diff_op,
            &Add::new(vec![
                AoTwoElecEnergy::builder(
                    TEST_OPER_NAME,
                    inner_density.differentiate(p.clone()).unwrap()
                )
                .outer_density(outer_density.clone())
                .dependencies(deps.clone())
                .build()
                .unwrap(),
                AoTwoElecEnergy::builder(TEST_OPER_NAME, inner_density.clone())
                    .outer_density(outer_density.differentiate(p).unwrap())
                    .dependencies(deps.clone())
                    .build()
                    .unwrap(),
            ])
            .unwrap()
        );
    }

    #[test]
    fn test_serialization() {
        let op =
            make_ao_two_elec_energy("", Some(make_wfn_parameter("")), Some(make_wfn_parameter("")));
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
        let op1 = make_ao_two_elec_energy(
            TEST_OPER_NAME,
            Some(inner_density.clone()),
            Some(outer_density.clone()),
        );

        assert!(is_expr_type::<AoTwoElecEnergy>(&op1));
        assert!(!is_zero_expr(&op1, None));
        assert!(!is_one_expr(&op1, None));

        let op2 = make_ao_two_elec_energy(
            TEST_OPER_NAME,
            Some(inner_density.clone()),
            Some(outer_density.clone()),
        );
        let op3 = make_ao_two_elec_energy(
            TEST_OPER_NAME,
            Some(outer_density.clone()),
            Some(inner_density.clone()),
        );
        let op4 =
            make_ao_two_elec_energy("", Some(inner_density.clone()), Some(outer_density.clone()));
        let op5 =
            make_ao_two_elec_energy(TEST_OPER_NAME, Some(make_wfn_parameter("density")), None);

        let op = downcast_from_arc::<AoTwoElecEnergy>(&op1).unwrap();
        let op6 = op
            .with_inner_density(outer_density.clone())
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
