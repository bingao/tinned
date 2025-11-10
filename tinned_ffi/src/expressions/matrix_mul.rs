use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::MatrixMul;
use tinned::public::generic_error;

use crate::c_support::{with_downcast_expr, with_downcast_val};
use crate::core::{
    ExprBox, ExprHandle, ExprSlice, TinnedErrorBox, expr_vec_from_slice, tinned_error_new,
};

#[ffi_export]
pub extern "C" fn tinned_matrix_mul_new(
    terms: ExprSlice<'_>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let terms_vec = match expr_vec_from_slice(terms, "tinned_matrix_mul_new") {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            return None;
        },
    };

    match MatrixMul::new(terms_vec) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

impl_expr_getters!(
    MatrixMul;
    tinned_matrix_mul_coefficient => |mul| Ok(Arc::clone(mul.coefficient())),
);

impl_val_getters!(
    MatrixMul;
    tinned_matrix_mul_factors_count: usize => |mul| mul.factors().len(); default = 0,
);

// Return the i-th factor (cloned). Caller must free the returned ExprBox.
impl_expr_index_getter!(tinned_matrix_mul_factor_at : MatrixMul => factors);
