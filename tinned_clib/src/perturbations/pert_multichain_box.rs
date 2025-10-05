use std::sync::Arc;

use tinned::perturbations::PertMultichain;

use crate::c_support::with_box_or_err;
use crate::core::TinnedErrorBox;

pub struct PertMultichainBox {
    inner: Arc<PertMultichain>,
}

impl PertMultichainBox {
    #[inline]
    pub(crate) fn arc_clone(&self) -> Arc<PertMultichain> {
        Arc::clone(&self.inner)
    }

    #[inline]
    pub(crate) fn as_ref(&self) -> &PertMultichain {
        self.inner.as_ref()
    }

    #[inline]
    pub(crate) fn as_arc_mut(&mut self) -> &mut Arc<PertMultichain> {
        &mut self.inner
    }
}

#[inline]
pub(crate) fn pert_multichain_box_from<T: Into<Arc<PertMultichain>>>(
    chain: T,
) -> *mut PertMultichainBox {
    Box::into_raw(Box::new(PertMultichainBox {
        inner: chain.into(),
    }))
}

// Borrows `&PertMultichain` or sets an error.
#[inline]
pub(crate) fn with_pert_multichain_or_err<R>(
    h: *const PertMultichainBox,
    out_err: *mut *mut TinnedErrorBox,
    caller: &'static str,
    f: impl FnOnce(&PertMultichain) -> R,
) -> Option<R> {
    with_box_or_err::<PertMultichainBox, PertMultichain, R>(
        h,
        out_err,
        caller,
        "PertMultichainBox",
        PertMultichainBox::as_ref,
        f,
    )
}
