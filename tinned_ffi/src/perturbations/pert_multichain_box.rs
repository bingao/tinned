use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::perturbations::PertMultichain;

/// An *opaque* handle that C can only pass around
#[derive_ReprC]
#[repr(opaque)]
pub struct PertMultichainHandle {
    inner: Arc<PertMultichain>,
}

/// Owned box by C after return
pub type PertMultichainBox = repr_c::Box<PertMultichainHandle>;

impl PertMultichainHandle {
    #[inline]
    pub(crate) fn new(p: Arc<PertMultichain>) -> Self {
        Self {
            inner: p,
        }
    }

    #[inline]
    pub(crate) fn as_ref(&self) -> &PertMultichain {
        &*self.inner
    }

    #[inline]
    pub(crate) fn clone_arc(&self) -> Arc<PertMultichain> {
        Arc::clone(&self.inner)
    }

    // Get a unique mutable reference to the chain if the Arc is uniquely owned. Used by `tinned_pert_multichain_insert`.
    #[inline]
    pub(crate) fn get_mut(&mut self) -> Option<&mut PertMultichain> {
        Arc::get_mut(&mut self.inner)
    }
}
