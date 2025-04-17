/// Represents exchange-correlation energy in a grid-based formulation.
/// xc_energy represents XC energy or its derivatives evaluated at grid points.
use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::utils::{build_xc_density, downcast_from_ref, intern_expr, validate_xc_inputs};

fn build_xc_energy(
    grid_weight: Arc<dyn Expr>,
    density_matrix: Arc<dyn Expr>,
    overlap_distribution: Arc<dyn Expr>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    let xc_density = build_xc_density("exc", density_matrix, overlap_distribution, 0)?;

    // - `Mul` for unperturbed or the first-order perturbed cases
    // - `Add` for higher-order perturbed case
    crate::expressions::Mul::new(vec![grid_weight, xc_density])
}

impl_exch_corr_type!(ExchCorrEnergy, ExchCorrEnergyBuilder, xc_energy, build_xc_energy);
impl_exch_corr_traits!(ExchCorrEnergy, xc_energy, true);

#[cfg(test)]
const DEFAULT_FUN_NAME: &str = "Exc[rho]";

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::non_elec_function::test_utils::make_non_elec_function;
    use crate::expressions::one_elec_operator::test_utils::make_one_elec_operator;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;

    #[inline]
    pub fn make_exch_corr_energy(
        name: impl Into<String>,
        grid_weight: Option<Arc<dyn Expr>>,
        density_matrix: Option<Arc<dyn Expr>>,
        overlap_distribution: Option<Arc<dyn Expr>>,
    ) -> Arc<dyn Expr> {
        let name: String = name.into();
        let weight = grid_weight.unwrap_or_else(|| make_non_elec_function(""));
        let dens = density_matrix.unwrap_or_else(|| make_wfn_parameter(""));
        let overlap = overlap_distribution.unwrap_or_else(|| make_one_elec_operator(""));
        if name.is_empty() {
            ExchCorrEnergy::builder(
                random_alphanumeric(DEFAULT_FUN_NAME.len() as u32 + 1),
                weight,
                dens,
                overlap,
            )
            .build()
            .unwrap()
        } else {
            ExchCorrEnergy::builder(name, weight, dens, overlap).build().unwrap()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;
    use crate::expressions::non_elec_function::test_utils::make_non_elec_function;
    use crate::expressions::one_elec_operator::test_utils::make_one_elec_operator;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_struct_safety!(ExchCorrEnergy);

    test_thread_interning!({
        make_exch_corr_energy(
            DEFAULT_FUN_NAME,
            Some(make_non_elec_function("weight")),
            Some(make_wfn_parameter("density")),
            Some(make_one_elec_operator("Omega")),
        )
    });

    #[test]
    fn test_impl_expr() {
        let weight = make_non_elec_function("");
        let density = make_wfn_parameter("");
        let overlap = make_one_elec_operator("");
        let op1 = make_exch_corr_energy(
            DEFAULT_FUN_NAME,
            Some(weight.clone()),
            Some(density.clone()),
            Some(overlap.clone()),
        );

        let op = downcast_from_arc::<ExchCorrEnergy>(&op1).unwrap();
        let xc_energy = build_xc_energy(weight.clone(), density.clone(), overlap.clone()).unwrap();
        assert_eq!(
            op,
            &ExchCorrEnergy {
                name: DEFAULT_FUN_NAME.into(),
                grid_weight: weight.clone(),
                density_matrix: density.clone(),
                overlap_distribution: overlap.clone(),
                xc_energy: xc_energy.clone(),
                derivative: PertMultichain::new(),
            }
        );

        assert_eq!(op.name(), DEFAULT_FUN_NAME);
        assert_eq!(op.grid_weight(), &weight);
        assert_eq!(op.density_matrix(), &density);
        assert_eq!(op.overlap_distribution(), &overlap);
        assert_eq!(op.xc_energy(), &xc_energy);

        assert_eq!(
            op1.hash_key(),
            format!(
                "ExchCorrEnergy({}; {}; {}; {}; [{}]; {})",
                DEFAULT_FUN_NAME,
                weight.hash_key(),
                density.hash_key(),
                overlap.hash_key(),
                PertMultichain::new().hash_key(),
                xc_energy.hash_key(),
            )
        );
        assert!(op1.is_scalar());
        assert_eq!(format!("{}", op1), format!("{}[{}]", DEFAULT_FUN_NAME, xc_energy));

        let op2 = make_exch_corr_energy(
            DEFAULT_FUN_NAME,
            Some(weight.clone()),
            Some(density.clone()),
            Some(overlap.clone()),
        );
        let op3 = make_exch_corr_energy(
            "",
            Some(weight.clone()),
            Some(density.clone()),
            Some(overlap.clone()),
        );
        let op4 = make_exch_corr_energy(
            DEFAULT_FUN_NAME,
            None,
            Some(density.clone()),
            Some(overlap.clone()),
        );
        let op5 =
            make_exch_corr_energy(DEFAULT_FUN_NAME, Some(weight.clone()), None, Some(overlap));
        let op6 = make_exch_corr_energy(DEFAULT_FUN_NAME, Some(weight), Some(density), None);

        assert_eq!(&op1, &op2);
        assert_ne!(&op1, &op3);
        assert_ne!(&op1, &op4);
        assert_ne!(&op1, &op5);
        assert_ne!(&op1, &op6);
    }

    #[test]
    fn test_serialization() {
        let op = make_exch_corr_energy("", None, None, None);
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let weight = make_non_elec_function("");
        let density = make_wfn_parameter("");
        let overlap = make_one_elec_operator("");
        let op1 = make_exch_corr_energy(
            DEFAULT_FUN_NAME,
            Some(weight.clone()),
            Some(density.clone()),
            Some(overlap.clone()),
        );
        let op2 = make_exch_corr_energy(
            DEFAULT_FUN_NAME,
            Some(weight.clone()),
            Some(density.clone()),
            Some(overlap.clone()),
        );
        let op3 = make_exch_corr_energy(
            "",
            Some(weight.clone()),
            Some(density.clone()),
            Some(overlap.clone()),
        );
        let op4 = make_exch_corr_energy(
            DEFAULT_FUN_NAME,
            None,
            Some(density.clone()),
            Some(overlap.clone()),
        );
        let op5 =
            make_exch_corr_energy(DEFAULT_FUN_NAME, Some(weight.clone()), None, Some(overlap));
        let op6 = make_exch_corr_energy(DEFAULT_FUN_NAME, Some(weight), Some(density), None);

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));
        assert!(!Arc::ptr_eq(&op1, &op4));
        assert!(!Arc::ptr_eq(&op1, &op5));
        assert!(!Arc::ptr_eq(&op1, &op6));

        assert!(is_expr_type::<ExchCorrEnergy>(&op1));
        assert!(!is_zero_expr(&op1));
        assert!(!is_one_expr(&op1));
    }
}
