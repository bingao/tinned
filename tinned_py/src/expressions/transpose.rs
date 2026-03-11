use pyo3::prelude::*;

use tinned::Transpose;

impl_unary_expr_interface!(
    expr_ty = Transpose,
    new_fn = transpose_new,
    new_doc = "Create a Transpose expression from an argument.",
    argument_fn = transpose_argument,
    argument_doc = "Return the argument of a Transpose expression.",
    register_fn = register
);
