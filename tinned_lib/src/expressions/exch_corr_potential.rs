/// Represents exchange-correlation potential in a grid-based formulation.
/// xc_potential represents XC potential or its derivatives evaluated at grid points.
use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::internal::build_xc_density;

fn build_xc_potential(
    grid_weight: Arc<dyn Expr>,
    density_matrix: Arc<dyn Expr>,
    overlap_distribution: Arc<dyn Expr>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    let xc_density = build_xc_density("vxc", density_matrix, overlap_distribution.clone(), 1)?;

    // - `MatrixMul` for unperturbed case, or when the generalized overlap
    //   distribution does not depend on applied perturbation(s)
    // - `MatrixAdd` for perturbed case in particular the generalized
    //   overlap distribution depends on applied perturbation(s)
    crate::expressions::MatrixMul::new(vec![
        crate::expressions::Mul::new(vec![grid_weight, xc_density])?,
        overlap_distribution,
    ])
}

impl_exch_corr_type!(
    ExchCorrPotential,
    ExchCorrPotentialBuilder,
    xc_potential,
    build_xc_potential,
    false
);
impl_exch_corr_traits!(ExchCorrPotential, xc_potential, false);

#[cfg(test)]
pub mod test_utils {
    use super::*;

    pub const TEST_FUNC_NAME: &str = "Vxc[rho]";

    impl_exch_corr_test_utils!(ExchCorrPotential, TEST_FUNC_NAME, make_exch_corr_potential);
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;

    test_exch_corr!(
        ExchCorrPotential,
        TEST_FUNC_NAME,
        make_exch_corr_potential,
        xc_potential,
        build_xc_potential,
        false
    );
}
