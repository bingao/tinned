use pyo3::prelude::*;

impl_addition_interface!(
    new_fn = matrix_add_new,
    new_fn_doc = "Create a MatrixAdd expression from a list of terms.",
    terms_fn = matrix_add_terms,
    terms_fn_doc = "Return the terms of a MatrixAdd expression as a list.",
    register_fn = register,
    expr_ty = tinned::MatrixAdd,
    downcast_err_msg = "matrix_add_terms() expected a MatrixAdd expression"
);
