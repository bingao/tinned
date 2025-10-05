use std::sync::Arc;

use tinned::perturbations::Perturbation;

use crate::c_support::{vec_arc_from_ptrs, with_box_or_err};
use crate::core::TinnedErrorBox;

pub struct PerturbationBox {
    inner: Arc<Perturbation>,
}

impl PerturbationBox {
    #[inline]
    pub(crate) fn arc_clone(&self) -> Arc<Perturbation> {
        Arc::clone(&self.inner)
    }

    #[inline]
    pub(crate) fn as_ref(&self) -> &Perturbation {
        self.inner.as_ref()
    }
}

// Allocates from an existing `Arc<Perturbation>`.
#[inline]
pub(crate) fn perturbation_box_from(pert: Arc<Perturbation>) -> *mut PerturbationBox {
    Box::into_raw(Box::new(PerturbationBox {
        inner: pert,
    }))
}

// Borrows `&Perturbation` or sets an error.
#[inline]
pub(crate) fn with_perturbation_or_err<R>(
    h: *const PerturbationBox,
    out_err: *mut *mut TinnedErrorBox,
    caller: &'static str,
    f: impl FnOnce(&Perturbation) -> R,
) -> Option<R> {
    with_box_or_err::<PerturbationBox, Perturbation, R>(
        h,
        out_err,
        caller,
        "PerturbationBox",
        PerturbationBox::as_ref,
        f,
    )
}

#[inline]
pub(crate) unsafe fn vec_pert_from_ptrs(
    ptrs: *const *const PerturbationBox,
    count: usize,
    caller: &'static str,
    out_err: *mut *mut TinnedErrorBox,
) -> Option<Vec<Arc<Perturbation>>> {
    unsafe {
        vec_arc_from_ptrs::<PerturbationBox, Perturbation>(
            ptrs,
            count,
            caller,
            "Perturbation",
            out_err,
            PerturbationBox::arc_clone,
        )
    }
}
