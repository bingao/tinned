/// Represents exchange-correlation energy in a grid-based formulation.
/// xc_energy represents XC energy or its derivatives evaluated at grid points.
use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::internal::build_xc_density;

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

impl_exch_corr_type!(ExchCorrEnergy, ExchCorrEnergyBuilder, xc_energy, build_xc_energy, true);
impl_exch_corr_traits!(ExchCorrEnergy, xc_energy, true);

#[cfg(test)]
pub mod test_utils {
    use super::*;

    pub const TEST_FUNC_NAME: &str = "Exc[rho]";

    impl_exch_corr_test_utils!(ExchCorrEnergy, TEST_FUNC_NAME, make_exch_corr_energy);
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;

    test_exch_corr!(
        ExchCorrEnergy,
        TEST_FUNC_NAME,
        make_exch_corr_energy,
        xc_energy,
        build_xc_energy,
        true
    );
}
