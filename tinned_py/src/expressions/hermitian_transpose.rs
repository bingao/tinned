use pyo3::prelude::*;

impl_unary_expr_interface!(
    new_fn = hermitian_transpose_new,
    new_fn_doc = "Create a HermitianTranspose expression from an argument.",
    arg_fn = hermitian_transpose_argument,
    arg_fn_doc = "Return the argument of a HermitianTranspose expression.",
    register_fn = register,
    expr_ty = tinned::HermitianTranspose,
    downcast_err_msg = "hermitian_transpose_argument() expected a HermitianTranspose expression"
);
