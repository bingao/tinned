use pyo3::prelude::*;

use tinned::NonElecFunction;

impl_nullary_expr_interface!(
    expr_ty = NonElecFunction,
    has_deps = true,
    new_fn = non_elec_function_new,
    name_fn = non_elec_function_name,
    is_perturbing_fn = non_elec_function_is_perturbing,
    derivative_fn = non_elec_function_derivative,
    deps_fn = non_elec_function_dependencies,
    indep_perts_fn = non_elec_function_independent_perturbations,
    register_fn = register
);
