use pyo3::prelude::*;

use tinned::TwoElecMatrix;

impl_nullary_expr_interface!(
    expr_ty = TwoElecMatrix,
    has_deps = true,
    new_fn = two_elec_matrix_new,
    name_fn = two_elec_matrix_name,
    is_perturbing_fn = two_elec_matrix_is_perturbing,
    derivative_fn = two_elec_matrix_derivative,
    deps_fn = two_elec_matrix_dependencies,
    indep_perts_fn = two_elec_matrix_independent_perturbations,
    register_fn = register
);
