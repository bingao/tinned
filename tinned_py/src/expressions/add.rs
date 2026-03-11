use pyo3::prelude::*;

use tinned::Add;

impl_addition_interface!(
    expr_ty = Add,
    new_fn = add_new,
    new_doc = "Create an Add expression from a list of terms.",
    terms_fn = add_terms,
    terms_doc = "Return the terms of an Add expression as a list.",
    register_fn = register
);
