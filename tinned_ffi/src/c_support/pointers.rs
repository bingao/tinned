use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::core::TinnedError;
use tinned::public::generic_error;

/// Validate an optional handle, then convert `&H` -> `T` (clone/produce).
/// Returns `None` and sets `out_err` if the handle is NULL.
#[inline]
pub fn try_from_handle<H, T>(
    h: Option<&H>,
    caller: &'static str,
    label: &'static str,
    to_value: impl FnOnce(&H) -> T,
) -> Result<T, TinnedError> {
    let h = h.ok_or_else(|| generic_error(format!("Null {label} pointer in {caller}"), None))?;
    Ok(to_value(h))
}

/// Validate an optional handle, map `&H` -> `&T`, then run `f(&T) -> R`.
/// Returns `None` and sets `out_err` if the handle is NULL.
#[inline]
pub fn try_with_handle<H, R>(
    h: Option<&H>,
    caller: &'static str,
    label: &'static str,
    f: impl FnOnce(&H) -> Result<R, TinnedError>,
) -> Result<R, TinnedError> {
    let h = h.ok_or_else(|| generic_error(format!("Null {label} pointer in {caller}"), None))?;
    f(h)
}

// Converts an array of `H` to `Vec<Arc<T>>`.
#[inline]
pub fn try_from_slice<H, T: ?Sized>(
    slice: c_slice::Ref<'_, *const H>,
    caller: &'static str,
    label: &'static str,
    mut to_arc: impl FnMut(&H) -> Arc<T>,
) -> Result<Vec<Arc<T>>, TinnedError> {
    let raw = slice.as_ref();
    let mut out = Vec::with_capacity(raw.len());
    for (i, &hptr) in raw.iter().enumerate() {
        if hptr.is_null() {
            return Err(generic_error(format!("Null {label} at index {i} in {caller}"), None));
        }
        let h: &H = unsafe { &*hptr };
        out.push(to_arc(h));
    }
    Ok(out)
}
