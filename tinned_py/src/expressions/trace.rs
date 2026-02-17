use pyo3::prelude::*;

impl_unary_expr_interface!(
    new_fn = trace_new,
    new_fn_doc = "Create a Trace expression from an argument.",
    arg_fn = trace_argument,
    arg_fn_doc = "Return the argument of a Trace expression.",
    register_fn = register,
    expr_ty = tinned::Trace,
    downcast_err_msg = "trace_argument() expected a Trace expression"
);
