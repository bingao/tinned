use pyo3::prelude::*;

impl_unary_expr_interface!(
    new_fn = transpose_new,
    new_fn_doc = "Create a Transpose expression from an argument.",
    arg_fn = transpose_argument,
    arg_fn_doc = "Return the argument of a Transpose expression.",
    register_fn = register,
    expr_ty = tinned::Transpose,
    downcast_err_msg = "transpose_argument() expected a Transpose expression"
);
