/// Represents exchange-correlation energy in a grid-based formulation.
/// xc_energy represents XC energy or its derivatives evaluated at grid points.
use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::internal::{build_xc_density, intern_expr, validate_xc_inputs};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{downcast_from_ref, generic_expression_error};

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
const DEFAULT_FUNC_NAME: &str = "Exc[rho]";

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::non_elec_function::test_utils::make_non_elec_function;
    use crate::expressions::one_elec_operator::test_utils::make_one_elec_operator;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;

    impl_exch_corr_test_utils!(ExchCorrEnergy, DEFAULT_FUNC_NAME, make_exch_corr_energy);
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;
    use crate::expressions::non_elec_function::test_utils::make_non_elec_function;
    use crate::expressions::one_elec_operator::test_utils::make_one_elec_operator;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_exch_corr!(
        ExchCorrEnergy,
        DEFAULT_FUNC_NAME,
        make_exch_corr_energy,
        xc_energy,
        build_xc_energy,
        true,
    );
}
