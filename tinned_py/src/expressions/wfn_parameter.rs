use pyo3::prelude::*;

use tinned::WfnParameter;

impl_nullary_expr_interface!(
    expr_ty = WfnParameter,
    has_deps = false,
    new_fn = wfn_parameter_new,
    name_fn = wfn_parameter_name,
    derivative_fn = wfn_parameter_derivative,
    register_fn = register
);
