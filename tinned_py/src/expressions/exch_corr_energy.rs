use pyo3::prelude::*;

use tinned::ExchCorrEnergy;

impl_exch_corr_interface!(
    expr_ty = ExchCorrEnergy,
    build_fn = exch_corr_energy_build,
    name_fn = exch_corr_energy_name,
    grid_weight_fn = exch_corr_energy_grid_weight,
    density_matrix_fn = exch_corr_energy_density_matrix,
    overlap_distribution_fn = exch_corr_energy_overlap_distribution,
    grid_expr_fn = exch_corr_energy_xc_energy,
    derivative_fn = exch_corr_energy_derivative,
    register_fn = register,
    grid_expr_method = xc_energy
);
