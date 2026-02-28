use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{MatrixAdd, ResidueParameter, WfnParameter, ZeroOperator};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{
    differentiate_expr, downcast_from_arc, downcast_from_ref, expression_error,
    generic_expression_error, is_expr_type, unreachable_error,
};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
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
    pub fn with_density(&self, density: Arc<dyn Expr>) -> TwoElecOperatorBuilder {
        TwoElecOperatorBuilder {
            name: self.name.clone(),
            density,
            dependencies: self.dependencies.clone(),
            derivative: self.derivative.clone(),
        }
    }

    #[inline]
    fn with_derivative(&self, derivative: PertMultichain) -> TwoElecOperatorBuilder {
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
    pub fn derivative(mut self, derivative: PertMultichain) -> Self {
        self.derivative = derivative;
        self
    }

    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        if is_expr_type::<ZeroOperator>(&self.density) {
            return Ok(self.density);
        }

        if is_expr_type::<WfnParameter>(&self.density)
            || is_expr_type::<ResidueParameter>(&self.density)
        {
            if self.dependencies.is_subchain(&self.derivative) {
                Ok(crate::internal::intern_expr(Arc::new(TwoElecOperator {
                    name: self.name,
                    density: self.density,
                    dependencies: self.dependencies,
                    derivative: self.derivative,
                })))
            } else {
                Ok(ZeroOperator::new())
            }
        } else {
            Err(expression_error(
                "TwoElecOperatorBuilder::build() - density must be WfnParameter or ResidueParameter",
                &self.density,
                None,
            ))
        }
    }
}

impl ExprInternal for TwoElecOperator {
    impl_unary_expr_internal_methods!(
        TwoElecOperator,
        density,
        true,
        |this: &TwoElecOperator, arg| this.with_density(arg).build()
    );

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "TwoElecOperator({}; {}; [{}]; [{}])",
            self.name,
            self.density.hash_key(),
            self.dependencies.hash_key(),
            self.derivative.hash_key(),
        )
    }

    #[inline]
    fn total_order(&self) -> u32 {
        self.derivative.total_order()
    }

    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<TwoElecOperator>(other) {
            self.name == op.name
                && self.dependencies == op.dependencies
                && self.density.deep_eq_superchains(&op.density)
                && self.derivative.is_subchain(&op.derivative)
        } else {
            false
        }
    }

    #[inline]
    fn eq_by_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        // For unambiguous replacement, we require equality of density
        // matrices, and make replacement by considering only derivative of
        // electron repulsion integrals (ERIs).
        if let Some(op) = downcast_from_arc::<TwoElecOperator>(other) {
            self.name == op.name
                && &self.density == &op.density
                && self.dependencies == op.dependencies
                && self.derivative.is_subchain(&op.derivative)
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for TwoElecOperator {
    impl_unary_expr_common_methods!(
        TwoElecOperator,
        density,
        False,
        |this: &TwoElecOperator, arg| this.with_density(arg).build()
    );

    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_density = self.density.differentiate(s).map_err(|e| {
            generic_expression_error(
                "TwoElecOperator::differentiate() failed for density",
                self,
                Some(Box::new(e)),
            )
        })?;
        let term1 = self.with_density(diff_density).build()?;

        let new_deriv = self.derivative.with_added_perturbation(s);
        let term2 = self.with_derivative(new_deriv).build()?;

        if is_expr_type::<ZeroOperator>(&term2) {
            return Ok(term1);
        }

        MatrixAdd::new(vec![term1, term2])
    }
}

impl PartialEq for TwoElecOperator {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
            && &self.density == &other.density
            && self.dependencies == other.dependencies
            && self.derivative == other.derivative
    }
}

impl Eq for TwoElecOperator {}

impl std::fmt::Display for TwoElecOperator {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}^{}[{}]", self.name, self.derivative, self.density)
    }
}

#[cfg(test)]
const DEFAULT_OPER_NAME: &str = "op(2el)";

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::pert_multichain::test_utils::{
        make_pert_multichain, make_super_multichain,
    };

    #[inline]
    pub fn make_two_elec_operator(
        name: impl Into<String>,
        density: Option<Arc<dyn Expr>>,
    ) -> Arc<dyn Expr> {
        let name: String = name.into();
        let dens = density.unwrap_or_else(|| make_wfn_parameter(""));
        if name.is_empty() {
            let deriv = make_pert_multichain(2u32, 10u32, 1u32, 10u32);
            let deps = make_super_multichain(&deriv, 1u32);
            TwoElecOperator::builder(random_alphanumeric(DEFAULT_OPER_NAME.len() as u32 + 1), dens)
                .dependencies(deps)
                .derivative(deriv)
                .build()
                .unwrap()
        } else {
            let deriv = make_pert_multichain(0u32, 0u32, 1u32, 0u32);
            let deps = make_super_multichain(&deriv, 1u32);
            TwoElecOperator::builder(name, dens)
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
    use crate::public::{is_one_expr, is_zero_expr};

    test_struct_safety!(TwoElecOperator);

    test_thread_interning!({
        make_two_elec_operator(DEFAULT_OPER_NAME, Some(make_wfn_parameter("density")))
    });

    #[test]
    fn test_impl_expr() {
        let density = make_wfn_parameter("");
        let deriv = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
        let op0 = TwoElecOperator::builder(DEFAULT_OPER_NAME, density.clone())
            .derivative(deriv.clone())
            .build()
            .unwrap();

        assert!(is_zero_expr(&op0, None));

        let deps = make_super_multichain(&deriv, 1u32);
        let op1 = TwoElecOperator::builder(DEFAULT_OPER_NAME, density.clone())
            .dependencies(deps.clone())
            .derivative(deriv.clone())
            .build()
            .unwrap();

        let op = downcast_from_arc::<TwoElecOperator>(&op1).unwrap();
        assert_eq!(
            op,
            &TwoElecOperator {
                name: DEFAULT_OPER_NAME.into(),
                density: density.clone(),
                dependencies: deps.clone(),
                derivative: deriv.clone()
            }
        );

        assert_eq!(op.name(), DEFAULT_OPER_NAME);
        assert_eq!(op.density(), &density);
        assert_eq!(op.dependencies(), &deps);
        assert_eq!(op.derivative(), &deriv);

        let mut op2 = op.with_density(density.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        op2 = op.with_derivative(deriv.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        assert_eq!(
            op1.hash_key(),
            format!(
                "TwoElecOperator({}; {}; [{}]; [{}])",
                DEFAULT_OPER_NAME,
                density.hash_key(),
                deps.hash_key(),
                deriv.hash_key(),
            )
        );
        assert!(!op1.is_scalar());
        assert_eq!(format!("{}", op1), format!("{}^{}[{}]", DEFAULT_OPER_NAME, deriv, density));

        let op3 = TwoElecOperator::builder(DEFAULT_OPER_NAME, density.clone())
            .dependencies(deps.clone())
            .derivative(deriv.clone())
            .build()
            .unwrap();
        let op4 = TwoElecOperator::builder(
            random_alphanumeric(DEFAULT_OPER_NAME.len() as u32 + 1),
            density.clone(),
        )
        .dependencies(deps.clone())
        .derivative(deriv.clone())
        .build()
        .unwrap();
        let op5 = TwoElecOperator::builder(DEFAULT_OPER_NAME, density.clone())
            .dependencies(deps.clone())
            .build()
            .unwrap();
        let op6 = TwoElecOperator::builder(DEFAULT_OPER_NAME, density.clone())
            .dependencies(make_super_multichain(&deriv, 2u32))
            .derivative(deriv.clone())
            .build()
            .unwrap();
        let op7 = TwoElecOperator::builder(DEFAULT_OPER_NAME, make_wfn_parameter("density"))
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
        let density = make_wfn_parameter("");
        let len_pert_name: u32 = 2;
        let deps = make_pert_multichain(len_pert_name, 8u32, 1u32, 10u32);
        let op = TwoElecOperator::builder(DEFAULT_OPER_NAME, density.clone())
            .dependencies(deps.clone())
            .build()
            .unwrap();

        let mut p: Arc<Perturbation> = deps.keys().first().cloned().unwrap();
        let mut diff_op = op.differentiate(&p).unwrap();
        let mut deriv = PertMultichain::new();
        deriv.insert(&p);

        assert_eq!(
            &diff_op,
            &MatrixAdd::new(vec![
                TwoElecOperator::builder(DEFAULT_OPER_NAME, density.clone())
                    .dependencies(deps.clone())
                    .derivative(deriv.clone())
                    .build()
                    .unwrap(),
                TwoElecOperator::builder(DEFAULT_OPER_NAME, density.differentiate(&p).unwrap())
                    .dependencies(deps.clone())
                    .build()
                    .unwrap(),
            ])
            .unwrap()
        );

        p = make_perturbation_symbol(len_pert_name + 1u32, 4u32);
        diff_op = op.differentiate(&p).unwrap();

        assert_eq!(
            &diff_op,
            &TwoElecOperator::builder(DEFAULT_OPER_NAME, density.differentiate(&p).unwrap())
                .dependencies(deps.clone())
                .build()
                .unwrap()
        );
    }

    #[test]
    fn test_serialization() {
        let op = make_two_elec_operator("", None);
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let density = make_wfn_parameter("");
        let op1 = make_two_elec_operator(DEFAULT_OPER_NAME, Some(density.clone()));

        assert!(is_expr_type::<TwoElecOperator>(&op1));
        assert!(!is_zero_expr(&op1, None));
        assert!(!is_one_expr(&op1, None));

        let op2 = make_two_elec_operator(DEFAULT_OPER_NAME, Some(density));
        let op3 = make_two_elec_operator("", Some(make_wfn_parameter("density")));
        let op4 = make_two_elec_operator(DEFAULT_OPER_NAME, Some(make_wfn_parameter("density")));

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));
        assert!(!Arc::ptr_eq(&op1, &op4));
    }
}
