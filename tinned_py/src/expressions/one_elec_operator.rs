use pyo3::prelude::*;

use tinned::OneElecOperator;

impl_nullary_expr_interface!(
    expr_ty = OneElecOperator,
    has_deps = true,
    new_fn = one_elec_operator_new,
    name_fn = one_elec_operator_name,
    derivative_fn = one_elec_operator_derivative,
    deps_fn = one_elec_operator_dependencies,
    register_fn = register
);
