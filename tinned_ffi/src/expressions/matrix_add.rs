use safer_ffi::prelude::*;

use tinned::expressions::MatrixAdd;

use crate::c_support::ffi_map_expr_as_exprvec;
use crate::core::{
    ExprBox, ExprHandle, ExprSlice, TinnedErrorBox, expr_vec_from_slice, tinned_error_new,
};

#[ffi_export]
pub extern "C" fn tinned_matrix_add_new(
    terms: &ExprSlice,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let terms_vec = match expr_vec_from_slice(terms, "tinned_matrix_add_new") {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            return None;
        },
    };

    match MatrixAdd::new(terms_vec) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Returns a cloned vector of terms
impl_vec_getter!(
    MatrixAdd,
    tinned_matrix_add_terms,
    repr_c::Vec<ExprBox>,
    ffi_map_expr_as_exprvec,
    |add| add.terms()
);
