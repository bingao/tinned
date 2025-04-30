use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{MatrixMul, OneElecOperator, TemporumOperator, ZeroOperator};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::utils::{downcast_from_ref, generic_expression_error, intern_expr, is_expr_type};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct TemporumOverlap {
    braket: Arc<dyn Expr>,
    dependencies: PertMultichain,
    derivative: PertMultichain,
}

impl TemporumOverlap {
    // Note: `dependencies` is the perturbation dependencies of Sb and Sk
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

// Helper function to build `braket` with given dependencies
fn build_braket(deps: &PertMultichain) -> Result<Arc<dyn Expr>, TinnedError> {
    let bra = OneElecOperator::builder("Sb").dependencies(deps.clone()).build()?;
    let dt_bra = TemporumOperator::builder(bra).is_forward(false).build()?;

    let ket = OneElecOperator::builder("Sk").dependencies(deps.clone()).build()?;
    let dt_ket = TemporumOperator::builder(ket).is_forward(true).build()?;

    MatrixMul::new(vec![dt_bra, dt_ket])
}

#[derive(Debug)]
pub struct TemporumOverlapBuilder {
    dependencies: PertMultichain,
}

impl TemporumOverlapBuilder {
    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(intern_expr(Arc::new(TemporumOverlap {
            braket: build_braket(&self.dependencies)?,
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
            self.dependencies.hash_key(),
            self.derivative.hash_key(),
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
        let diff_braket = self.braket.differentiate(s).map_err(|e| {
            generic_expression_error("Differentiation failed", self, Some(Box::new(e)))
        })?;

        if is_expr_type::<ZeroOperator>(&diff_braket) {
            return Ok(diff_braket);
        }

        let new_deriv = self.derivative.clone_with_insert(s);

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
        write!(f, "op(T)^{}", self.derivative)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perturbations::pert_multichain::test_utils::{
        make_pert_multichain, make_super_multichain,
    };
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::utils::{downcast_from_arc, is_one_expr, is_zero_expr};

    test_struct_safety!(TemporumOverlap);

    test_thread_interning!(
        TemporumOverlap::builder(make_pert_multichain(0u32, 0u32, 1u32, 0u32)).build().unwrap()
    );

    #[test]
    fn test_impl_expr() {
        let deps = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
        let op1 = TemporumOverlap::builder(deps.clone()).build().unwrap();

        let op = downcast_from_arc::<TemporumOverlap>(&op1).unwrap();
        assert_eq!(
            op,
            &TemporumOverlap {
                braket: build_braket(&deps).unwrap(),
                dependencies: deps.clone(),
                derivative: PertMultichain::new(),
            }
        );

        assert_eq!(op.dependencies(), &deps);

        assert_eq!(
            op1.hash_key(),
            format!(
                "TemporumOverlap([{}]; [{}])",
                deps.hash_key(),
                PertMultichain::new().hash_key(),
            )
        );
        assert!(!op1.is_scalar());
        assert_eq!(format!("{}", op1), format!("op(T)^{}", PertMultichain::new()));

        let op2 = TemporumOverlap::builder(make_super_multichain(&deps, 1u32)).build().unwrap();

        assert_ne!(&op1, &op2);
    }

    #[test]
    fn test_differentiation() {
        let len_pert_name: u32 = 2;
        let deps = make_pert_multichain(len_pert_name, 8u32, 1u32, 10u32);
        let op = TemporumOverlap::builder(deps.clone()).build().unwrap();

        let p: Arc<Perturbation> = deps.keys().first().cloned().unwrap();
        let mut diff_op = op.differentiate(&p).unwrap();
        let mut deriv = PertMultichain::new();
        deriv.insert(&p);

        let diff_cast = downcast_from_arc::<TemporumOverlap>(&diff_op).unwrap();

        assert_eq!(diff_cast.derivative(), &deriv);

        let max_order = deps.get_order(&p);
        for _ in 1..=2 * max_order + 1 {
            diff_op = diff_op.differentiate(&p).unwrap();
        }

        assert_eq!(&diff_op, &ZeroOperator::new());

        assert!(is_zero_expr(
            &op.differentiate(&make_perturbation_symbol(len_pert_name + 1u32, 4u32)).unwrap(),
            None
        ));
    }

    #[test]
    fn test_serialization() {
        let op = TemporumOverlap::builder(make_pert_multichain(2u32, 8u32, 1u32, 10u32))
            .build()
            .unwrap();
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let deps = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
        let op1 = TemporumOverlap::builder(deps.clone()).build().unwrap();

        assert!(is_expr_type::<TemporumOverlap>(&op1));
        assert!(!is_zero_expr(&op1, None));
        assert!(!is_one_expr(&op1, None));

        let op2 = TemporumOverlap::builder(deps.clone()).build().unwrap();
        let op3 = TemporumOverlap::builder(make_super_multichain(&deps, 1u32)).build().unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));
    }
}
