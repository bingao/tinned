use std::sync::Arc;

use tinned::public::generic_error;

use crate::core::{TinnedErrorBox, set_out_err};

// Borrows `&T` safely from a raw handle. Returns `None` if `h` is `NULL`.
#[inline]
pub(crate) fn with_box_or_err<H, T: ?Sized, R>(
    h: *const H,
    out_err: *mut *mut TinnedErrorBox,
    caller: &'static str,
    label: &'static str,
    to_target: fn(&H) -> &T,
    f: impl FnOnce(&T) -> R,
) -> Option<R> {
    if h.is_null() {
        set_out_err(out_err, generic_error(format!("Null {label} pointer in {caller}"), None));
        return None;
    }

    let b: &H = unsafe { &*h };
    let target: &T = to_target(b);
    Some(f(target))
}

// Converts an array of `H` to `Vec<Arc<T>>`.
// - `ptrs` must point to an array of `count` elements of type `*const H`,
//   properly aligned and alive for the duration of this call.
// - Each element may be null (we check and error out), otherwise must point to a valid `H`.
pub(crate) unsafe fn vec_arc_from_ptrs<H, T: ?Sized>(
    ptrs: *const *const H,
    count: usize,
    caller: &'static str,
    label: &'static str,
    out_err: *mut *mut TinnedErrorBox,
    to_arc: fn(&H) -> Arc<T>,
) -> Option<Vec<Arc<T>>> {
    if ptrs.is_null() {
        if count == 0 {
            return Some(Vec::new());
        }
        set_out_err(out_err, generic_error(format!("Null pointer array in {caller}"), None));
        return None;
    }

    let slice: &[*const H] = unsafe { std::slice::from_raw_parts(ptrs, count) };

    let mut out = Vec::with_capacity(count);
    for (i, &h) in slice.iter().enumerate() {
        if h.is_null() {
            set_out_err(
                out_err,
                generic_error(format!("Null {label} at index {i} in {caller}"), None),
            );
            return None;
        }
        // Borrows the handle
        let b: &H = unsafe { &*h };
        // `&H` -> `Arc<T>`
        let arc = to_arc(b);
        out.push(arc);
    }
    Some(out)
}
