/// Represents exchange-correlation potential in a grid-based formulation.
/// xc_potential represents XC potential or its derivatives evaluated at grid points.
use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::utils::{build_xc_density, downcast_from_ref, intern_expr, validate_xc_inputs};

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

impl_exch_corr_type!(ExchCorrPotential, ExchCorrPotentialBuilder, xc_potential, build_xc_potential);
impl_exch_corr_traits!(ExchCorrPotential, xc_potential, false);

#[cfg(test)]
const DEFAULT_FUNC_NAME: &str = "Vxc[rho]";

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::non_elec_function::test_utils::make_non_elec_function;
    use crate::expressions::one_elec_operator::test_utils::make_one_elec_operator;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;

    impl_exch_corr_test_utils!(ExchCorrPotential, DEFAULT_FUNC_NAME, make_exch_corr_potential);
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;
    use crate::expressions::non_elec_function::test_utils::make_non_elec_function;
    use crate::expressions::one_elec_operator::test_utils::make_one_elec_operator;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_exch_corr!(
        ExchCorrPotential,
        DEFAULT_FUNC_NAME,
        make_exch_corr_potential,
        xc_potential,
        build_xc_potential,
        false
    );
}
