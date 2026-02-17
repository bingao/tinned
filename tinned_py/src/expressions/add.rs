use pyo3::prelude::*;

impl_addition_interface!(
    new_fn = add_new,
    new_fn_doc = "Create an Add expression from a list of terms.",
    terms_fn = add_terms,
    terms_fn_doc = "Return the terms of an Add expression as a list.",
    register_fn = register,
    expr_ty = tinned::Add,
    downcast_err_msg = "add_terms() expected an Add expression"
);
