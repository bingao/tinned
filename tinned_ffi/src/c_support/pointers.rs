use safer_ffi::prelude::*;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
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
pub fn try_vec_from_slice<H, T: ?Sized>(
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

// Converts an array of `H` to `HashSet<Arc<T>>`.
#[inline]
pub fn try_set_from_slice<H, T: ?Sized + Eq + Hash>(
    slice: c_slice::Ref<'_, *const H>,
    caller: &'static str,
    label: &'static str,
    mut to_arc: impl FnMut(&H) -> Arc<T>,
) -> Result<HashSet<Arc<T>>, TinnedError> {
    let raw = slice.as_ref();
    let mut out: HashSet<Arc<T>> = HashSet::with_capacity(raw.len());
    for (i, &hptr) in raw.iter().enumerate() {
        if hptr.is_null() {
            return Err(generic_error(format!("Null {label} at index {i} in {caller}"), None));
        }
        let h: &H = unsafe { &*hptr };
        out.insert(to_arc(h));
    }
    Ok(out)
}

// Converts two parallel arrays into `HashMap<Arc<K>, Arc<V>>`.
#[inline]
pub fn try_map_from_slices<HK, HV, K: ?Sized + Eq + Hash, V: ?Sized>(
    keys: c_slice::Ref<'_, *const HK>,
    values: c_slice::Ref<'_, *const HV>,
    caller: &'static str,
    key_label: &'static str,
    val_label: &'static str,
    mut to_key: impl FnMut(&HK) -> Arc<K>,
    mut to_val: impl FnMut(&HV) -> Arc<V>,
) -> Result<HashMap<Arc<K>, Arc<V>>, TinnedError> {
    let kraw = keys.as_ref();
    let vraw = values.as_ref();
    if kraw.len() != vraw.len() {
        return Err(generic_error(
            format!("Mismatched key/value lengths ({} vs {}) in {caller}", kraw.len(), vraw.len()),
            None,
        ));
    }
    let mut out: HashMap<Arc<K>, Arc<V>> = HashMap::with_capacity(kraw.len());

    for i in 0..kraw.len() {
        let kptr = kraw[i];
        let vptr = vraw[i];
        if kptr.is_null() {
            return Err(generic_error(format!("Null {key_label} at index {i} in {caller}"), None));
        }
        if vptr.is_null() {
            return Err(generic_error(format!("Null {val_label} at index {i} in {caller}"), None));
        }
        let kh: &HK = unsafe { &*kptr };
        let vh: &HV = unsafe { &*vptr };
        out.insert(to_key(kh), to_val(vh));
    }
    Ok(out)
}
