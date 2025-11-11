use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::core::TinnedError;
use tinned::perturbations::Perturbation;
use tinned::public::generic_error;

use crate::c_support::{
    tinned_string_from_cstr, tinned_string_to_cstr, try_from_handle, try_vec_from_slice,
    try_with_handle,
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

/// Borrowed slice of handles
pub type PerturbationSlice<'a> = c_slice::Ref<'a, *const PerturbationHandle>;

/// Turn a `PerturbationSlice` into `Vec<Arc<Perturbation>>`.
#[inline]
pub fn perturbation_vec_from_slice(
    slice: PerturbationSlice<'_>,
    caller: &'static str,
) -> Result<Vec<Arc<Perturbation>>, TinnedError> {
    // Reuse the same safety/validation logic as Expr via `try_vec_from_slice`
    try_vec_from_slice(slice, caller, "PerturbationHandle", |h: &PerturbationHandle| h.clone_arc())
}

// Free a perturbation (NULL-safe).
#[ffi_export]
pub fn tinned_perturbation_free(pert: Option<PerturbationBox>) {
    drop(pert);
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

// Get `name`; caller frees with tinned_string_free. NULL on error.
#[ffi_export]
pub extern "C" fn tinned_perturbation_name(
    h: Option<&PerturbationHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    match try_with_handle(h, "tinned_perturbation_name", "PerturbationHandle", |ph| {
        let p = ph.as_ref();
        Ok(tinned_string_to_cstr(p.name().to_string()))
    }) {
        Ok(s) => Some(s),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Get `frequency` (cloned). Caller must free the returned ExprBox. NULL on error.
#[ffi_export]
pub extern "C" fn tinned_perturbation_frequency(
    h: Option<&PerturbationHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    match try_with_handle(h, "tinned_perturbation_frequency", "PerturbationHandle", |ph| {
        let p = ph.as_ref();
        let freq = Arc::clone(p.frequency());
        Ok(ExprBox::new(ExprHandle::new(freq)))
    }) {
        Ok(b) => Some(b),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Display text; caller frees with tinned_string_free. NULL on error.
#[ffi_export]
pub extern "C" fn tinned_perturbation_display(
    h: Option<&PerturbationHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    match try_with_handle(h, "tinned_perturbation_display", "PerturbationHandle", |ph| {
        let p = ph.as_ref();
        Ok(tinned_string_to_cstr(format!("{}", p)))
    }) {
        Ok(s) => Some(s),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}
