use pyo3::prelude::*;

impl_nullary_expr_interface!(
    type_name = tinned::LagMultiplier,
    type_label = "LagMultiplier",
    has_deps = false,
    new_fn = lag_multiplier_new,
    name_fn = lag_multiplier_name,
    derivative_fn = lag_multiplier_derivative,
    register_fn = register
);
