use std::ptr::null_mut;

use tinned::expressions::Add;
use tinned::public::generic_error;

use crate::c_support::{with_downcast_expr_res, with_downcast_val};
use crate::core::{ExprBox, TinnedErrorBox, set_out_err, vec_expr_from_ptrs};

#[unsafe(no_mangle)]
pub extern "C" fn tinned_add_new(
    term_ptrs: *const *const ExprBox,
    term_count: usize,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
    let terms = unsafe {
        match vec_expr_from_ptrs(term_ptrs, term_count, "tinned_add_new", out_err) {
            Some(v) => v,
            None => return null_mut(),
        }
    };

    match Add::new(terms) {
        Ok(expr) => ExprBox::new(expr).into_raw(),
        Err(e) => {
            set_out_err(out_err, e);
            null_mut()
        },
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_add_terms_count(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> usize {
    with_downcast_val::<Add, usize>(h, out_err, "tinned_add_terms_count", |a| a.terms().len())
}

// Return the i-th term (cloned). Caller must unref.
#[unsafe(no_mangle)]
pub extern "C" fn tinned_add_term_at(
    h: *const ExprBox,
    i: usize,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
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
