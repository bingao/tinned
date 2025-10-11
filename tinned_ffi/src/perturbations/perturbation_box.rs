use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::core::TinnedError;
use tinned::perturbations::Perturbation;

use crate::c_support::try_from_slice;

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
    pub(crate) fn new(p: Arc<Perturbation>) -> Self {
        Self {
            inner: p,
        }
    }

    #[inline]
    pub(crate) fn as_ref(&self) -> &Perturbation {
        &*self.inner
    }

    #[inline]
    pub(crate) fn clone_arc(&self) -> Arc<Perturbation> {
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
    // Reuse the same safety/validation logic as Expr via `try_from_slice`
    try_from_slice(slice, caller, "PerturbationHandle", |h: &PerturbationHandle| h.clone_arc())
}
