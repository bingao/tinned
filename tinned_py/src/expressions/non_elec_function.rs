use pyo3::prelude::*;

impl_nullary_expr_interface!(
    type_name = tinned::NonElecFunction,
    type_label = "NonElecFunction",
    has_deps = true,
    new_fn = non_elec_function_new,
    name_fn = non_elec_function_name,
    derivative_fn = non_elec_function_derivative,
    deps_fn = non_elec_function_dependencies,
    register_fn = register
);
