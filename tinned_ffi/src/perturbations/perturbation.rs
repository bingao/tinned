use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::core::TinnedError;
use tinned::perturbations::Perturbation;
use tinned::public::generic_error;

use crate::c_support::{
    tinned_string_from_cstr, tinned_string_to_cstr, try_from_handle, try_set_from_slice,
    try_vec_from_slice, try_with_handle,
};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};

/// An *opaque* handle that C can only pass around
#[derive_ReprC]
#[repr(opaque)]
pub struct PerturbationHandle {
    inner: Arc<Perturbation>,
}

/// Owned box by C after return
pub type PerturbationBox = repr_c::Box<PerturbationHandle>;

impl PerturbationHandle {
    #[inline]
    pub fn new(p: Arc<Perturbation>) -> Self {
        Self {
            inner: p,
        }
    }

    #[inline]
    pub fn as_ref(&self) -> &Perturbation {
        &*self.inner
    }

    #[inline]
    pub fn clone_arc(&self) -> Arc<Perturbation> {
        Arc::clone(&self.inner)
    }
}

/// Free a perturbation (NULL-safe).
#[ffi_export]
pub fn tinned_perturbation_free(pert: Option<PerturbationBox>) {
    drop(pert);
}

/// Frees a vector of perturbation boxes returned by Rust.
///
/// This also drops all contained `PerturbationBox` elements.
/// The caller must not free the elements separately afterward.
#[ffi_export]
pub fn tinned_perturbation_vec_free(_v: repr_c::Vec<PerturbationBox>) {
    // Intentionally empty. Taking by value and returning lets _v drop here.
}

#[inline]
fn with_perturbation_ref<R>(
    h: Option<&PerturbationHandle>,
    caller: &'static str,
    f: impl FnOnce(&Perturbation) -> Result<R, TinnedError>,
) -> Result<R, TinnedError> {
    try_with_handle(h, caller, "PerturbationHandle", |ph| {
        let pert = ph.as_ref();
        f(pert)
    })
}

#[inline]
fn ffi_perturbation_return_val<R>(
    h: Option<&PerturbationHandle>,
    caller: &'static str,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    f: impl FnOnce(&Perturbation) -> Result<R, TinnedError>,
) -> R
where
    R: Default,
{
    match with_perturbation_ref(h, caller, f) {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            R::default()
        },
    }
}

// Create a new Perturbation.
#[ffi_export]
pub extern "C" fn tinned_perturbation_new(
    name: Option<char_p::Ref<'_>>,
    frequency: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<PerturbationBox> {
    let Some(name) = tinned_string_from_cstr(name) else {
        tinned_error_new(
            out_err,
            generic_error("Null perturbation name passed to tinned_perturbation_new", None),
        );
        return None;
    };

    let Some(freq) = frequency else {
        tinned_error_new(
            out_err,
            generic_error("Null frequency passed to tinned_perturbation_new", None),
        );
        return None;
    };

    let freq_arc = freq.clone_arc();
    let pert = Perturbation::new(name, freq_arc);
    Some(PerturbationBox::new(PerturbationHandle::new(pert)))
}

// Clone a perturbation (like Arc clone). Returns NULL on error / NULL input.
#[ffi_export]
pub fn tinned_perturbation_clone(
    h: Option<&PerturbationHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<PerturbationBox> {
    match try_from_handle(h, "tinned_perturbation_clone", "PerturbationHandle", |ph| {
        PerturbationBox::new(PerturbationHandle::new(ph.clone_arc()))
    }) {
        Ok(b) => Some(b),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Get `name`; caller frees with tinned_string_free. NULL on error.
#[ffi_export]
pub extern "C" fn tinned_perturbation_name(
    h: Option<&PerturbationHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    ffi_perturbation_return_val(h, "tinned_perturbation_name", out_err, |pert| {
        Ok(Some(tinned_string_to_cstr(pert.name())))
    })
}

// Get `frequency` (cloned). Caller must free the returned ExprBox. NULL on error.
#[ffi_export]
pub extern "C" fn tinned_perturbation_frequency(
    h: Option<&PerturbationHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    ffi_perturbation_return_val(h, "tinned_perturbation_frequency", out_err, |pert| {
        let freq = Arc::clone(pert.frequency());
        Ok(Some(ExprBox::new(ExprHandle::new(freq))))
    })
}

// Display text; caller frees with tinned_string_free. NULL on error.
#[ffi_export]
pub extern "C" fn tinned_perturbation_display(
    h: Option<&PerturbationHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    ffi_perturbation_return_val(h, "tinned_perturbation_display", out_err, |pert| {
        Ok(Some(tinned_string_to_cstr(format!("{}", pert))))
    })
}

// One (perturbation, max_order) entry
#[repr(C)]
#[derive_ReprC]
pub struct PerturbationEntry {
    perturbation: repr_c::Box<PerturbationHandle>,
    max_order: u32,
}

impl PerturbationEntry {
    #[inline]
    pub fn new(perturbation: repr_c::Box<PerturbationHandle>, max_order: u32) -> Self {
        Self {
            perturbation,
            max_order,
        }
    }

    #[inline]
    pub fn perturbation(&self) -> &repr_c::Box<PerturbationHandle> {
        &self.perturbation
    }

    #[inline]
    pub fn max_order(&self) -> u32 {
        self.max_order
    }
}

// FFI constructor for C
#[ffi_export]
pub extern "C" fn tinned_perturbation_entry_new(
    perturbation: repr_c::Box<PerturbationHandle>,
    max_order: u32,
) -> PerturbationEntry {
    PerturbationEntry::new(perturbation, max_order)
}

/// Borrowed slice of handles
#[repr(C)]
#[derive_ReprC]
pub struct PerturbationSlice {
    pub ptr: *const *const PerturbationHandle,
    pub len: usize,
}

/// Turn a `PerturbationSlice` into `Vec<Arc<Perturbation>>`.
#[inline]
pub fn perturbation_vec_from_slice(
    slice: &PerturbationSlice,
    caller: &'static str,
) -> Result<Vec<Arc<Perturbation>>, TinnedError> {
    try_vec_from_slice(
        slice.ptr,
        slice.len,
        caller,
        "PerturbationHandle",
        |h: &PerturbationHandle| h.clone_arc(),
    )
}

// Build a `HashSet<Arc<Perturbation>>` or `BTreeSet<Arc<Perturbation>>` from a `PerturbationSlice`.
#[inline]
pub fn perturbation_set_from_slice<S>(
    slice: &PerturbationSlice,
    caller: &'static str,
) -> Result<S, TinnedError>
where
    S: Default + Extend<Arc<Perturbation>>,
{
    try_set_from_slice(
        slice.ptr,
        slice.len,
        caller,
        "PerturbationHandle",
        |h: &PerturbationHandle| h.clone_arc(),
    )
}

/// Slice of `PerturbationEntry`
#[repr(C)]
#[derive_ReprC]
pub struct PerturbationEntrySlice {
    pub ptr: *const PerturbationEntry,
    pub len: usize,
}
