use pyo3::prelude::*;

use tinned::MatrixAdd;

impl_addition_interface!(
    expr_ty = MatrixAdd,
    new_fn = matrix_add_new,
    new_doc = "Create a MatrixAdd expression from a list of terms.",
    terms_fn = matrix_add_terms,
    terms_doc = "Return the terms of a MatrixAdd expression as a list.",
    register_fn = register
);
