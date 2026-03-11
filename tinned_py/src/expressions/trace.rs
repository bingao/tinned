use pyo3::prelude::*;

use tinned::Trace;

impl_unary_expr_interface!(
    expr_ty = Trace,
    new_fn = trace_new,
    new_doc = "Create a Trace expression from an argument.",
    argument_fn = trace_argument,
    argument_doc = "Return the argument of a Trace expression.",
    register_fn = register
);
