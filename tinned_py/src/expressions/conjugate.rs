use pyo3::prelude::*;

impl_unary_expr_interface!(
    new_fn = conjugate_new,
    new_fn_doc = "Create a Conjugate expression from an argument.",
    arg_fn = conjugate_argument,
    arg_fn_doc = "Return the argument of a Conjugate expression.",
    register_fn = register,
    expr_ty = tinned::Conjugate,
    downcast_err_msg = "conjugate_argument() expected a Conjugate expression"
);
