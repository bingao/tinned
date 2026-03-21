use pyo3::prelude::*;

use tinned::OneElecMatrix;

impl_nullary_expr_interface!(
    expr_ty = OneElecMatrix,
    has_deps = true,
    new_fn = one_elec_matrix_new,
    name_fn = one_elec_matrix_name,
    is_perturbing_fn = one_elec_matrix_is_perturbing,
    derivative_fn = one_elec_matrix_derivative,
    deps_fn = one_elec_matrix_dependencies,
    indep_perts_fn = one_elec_matrix_independent_perturbations,
    register_fn = register
);
