use safer_ffi::prelude::*;

use tinned::expressions::MatrixAdd;
use tinned::public::generic_error;

use crate::c_support::{with_downcast_expr, with_downcast_val};
use crate::core::{
    ExprBox, ExprHandle, ExprSlice, TinnedErrorBox, expr_vec_from_slice, tinned_error_new,
};

#[ffi_export]
pub extern "C" fn tinned_matrix_add_new(
    terms: ExprSlice<'_>,
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

impl_val_getters!(
    MatrixAdd;
    tinned_matrix_add_terms_count: usize => |add| add.terms().len(); default = 0,
);

// Return the i-th term (cloned). Caller must free the returned ExprBox.
impl_expr_index_getter!(tinned_matrix_add_term_at : MatrixAdd => terms);
