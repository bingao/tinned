use std::{os::raw::c_char, ptr::null_mut, sync::Arc};

use tinned::core::Expr;
use tinned::public::generic_error;

use crate::c_support::{cstr_to_string, to_cstring};
use crate::core::{ExprBox, TinnedErrorBox, expr_box_from, set_out_err, with_expr_or_err};

// Clones the `Arc` inside `ExprBox`, increase the strong count, and return a
// new pointer to a freshly boxed ExprBox that points to the same underlying
// Rust object. Must eventually call `tinned_expr_unref()` on the returned
// pointer.
#[unsafe(no_mangle)]
pub extern "C" fn tinned_expr_ref(h: *mut ExprBox) -> *mut ExprBox {
    if h.is_null() {
        return null_mut();
    }
    let expr = unsafe { &*h }.arc_clone();
    expr_box_from(expr)
}

// Drops `ExprBox`. That in turn drops one `Arc<dyn Expr>` strong reference.
// When the last reference goes away, the Rust object is deallocated. Call this
// exactly once for every owned pointer.
#[unsafe(no_mangle)]
pub extern "C" fn tinned_expr_unref(h: *mut ExprBox) {
    if h.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(h));
    }
}

// Public APIs for accessing functions of `Expr`

#[unsafe(no_mangle)]
pub extern "C" fn tinned_expr_is_scalar(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> bool {
    with_expr_or_err(h, out_err, "tinned_expr_is_scalar", |e| e.is_scalar()).unwrap_or(false)
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_expr_hash_key(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut c_char {
    with_expr_or_err(h, out_err, "tinned_expr_hash_key", |e| e.hash_key())
        .map(to_cstring)
        .unwrap_or(null_mut())
}

// Get printable expression text. C must free with `tinned_free_cstring()`.
#[unsafe(no_mangle)]
pub extern "C" fn tinned_expr_display(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut c_char {
    with_expr_or_err(h, out_err, "tinned_expr_display", |e| format!("{}", e))
        .map(to_cstring)
        .unwrap_or(null_mut())
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_expr_serialize_json(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut c_char {
    with_expr_or_err(h, out_err, "tinned_expr_serialize_json", |e| match serde_json::to_string(e) {
        Ok(s) => to_cstring(s),
        Err(err) => {
            set_out_err(
                out_err,
                generic_error("Failed to serialize expression to JSON", Some(Box::new(err))),
            );
            null_mut()
        },
    })
    .unwrap_or(null_mut())
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_expr_deserialize_json(json: *const c_char) -> *mut ExprBox {
    let s = match cstr_to_string(json) {
        Some(s) => s,
        None => return null_mut(),
    };
    match serde_json::from_str::<Arc<dyn Expr>>(&s) {
        Ok(arc) => expr_box_from(arc),
        Err(_) => null_mut(),
    }
}
