use pyo3::prelude::*;

use tinned::HermitianTranspose;

impl_unary_expr_interface!(
    expr_ty = HermitianTranspose,
    new_fn = hermitian_transpose_new,
    new_doc = "Create a HermitianTranspose expression from an argument.",
    argument_fn = hermitian_transpose_argument,
    argument_doc = "Return the argument of a HermitianTranspose expression.",
    register_fn = register
);
