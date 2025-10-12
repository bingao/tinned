use safer_ffi::prelude::*;

use tinned::expressions::Add;
use tinned::public::generic_error;

use crate::c_support::{with_downcast_expr_res, with_downcast_val};
use crate::core::{
    ExprBox, ExprHandle, ExprSlice, TinnedErrorBox, expr_vec_from_slice, tinned_error_new,
};

#[ffi_export]
pub extern "C" fn tinned_add_new(
    terms: ExprSlice<'_>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let terms_vec = match expr_vec_from_slice(terms, "tinned_add_new") {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            return None;
        },
    };

    match Add::new(terms_vec) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

#[ffi_export]
pub extern "C" fn tinned_add_terms_count(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> usize {
    with_downcast_val::<Add, usize>(h, out_err, "tinned_add_terms_count", |a| a.terms().len())
        .unwrap_or(0)
}

// Return the i-th term (cloned). Caller must free the returned ExprBox.
#[ffi_export]
pub extern "C" fn tinned_add_term_at(
    h: Option<&ExprHandle>,
    i: usize,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    with_downcast_expr_res::<Add>(h, out_err, "tinned_add_term_at", |a| {
        a.terms().get(i).cloned().ok_or_else(|| {
            generic_error(
                format!(
                    "Index {i} out of bounds (len = {}) in tinned_add_term_at",
                    a.terms().len()
                ),
                None,
            )
        })
    })
}
