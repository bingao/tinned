use pyo3::prelude::*;

use tinned::ExchCorrPotential;

impl_exch_corr_interface!(
    expr_ty = ExchCorrPotential,
    build_fn = exch_corr_potential_build,
    name_fn = exch_corr_potential_name,
    grid_weight_fn = exch_corr_potential_grid_weight,
    density_matrix_fn = exch_corr_potential_density_matrix,
    overlap_distribution_fn = exch_corr_potential_overlap_distribution,
    grid_expr_fn = exch_corr_potential_xc_potential,
    derivative_fn = exch_corr_potential_derivative,
    register_fn = register,
    grid_expr_method = xc_potential
);
