use pyo3::prelude::*;

use tinned::LagMultiplier;

impl_nullary_expr_interface!(
    expr_ty = LagMultiplier,
    has_deps = false,
    new_fn = lag_multiplier_new,
    name_fn = lag_multiplier_name,
    derivative_fn = lag_multiplier_derivative,
    register_fn = register
);
