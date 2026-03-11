use pyo3::prelude::*;

use tinned::Conjugate;

impl_unary_expr_interface!(
    expr_ty = Conjugate,
    new_fn = conjugate_new,
    new_doc = "Create a Conjugate expression from an argument.",
    argument_fn = conjugate_argument,
    argument_doc = "Return the argument of a Conjugate expression.",
    register_fn = register
);
