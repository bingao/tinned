use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::perturbations::PertMultichain;
use tinned::public::generic_error;

use crate::c_support::{tinned_string_to_cstr, try_from_handle, try_with_handle};
use crate::core::{TinnedErrorBox, tinned_error_new};
use crate::perturbations::{PerturbationHandle, PerturbationSlice, perturbation_vec_from_slice};

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
    pub fn new(p: Arc<PertMultichain>) -> Self {
        Self {
            inner: p,
        }
    }

    #[inline]
    pub fn as_ref(&self) -> &PertMultichain {
        &*self.inner
    }

    #[inline]
    pub fn clone_arc(&self) -> Arc<PertMultichain> {
        Arc::clone(&self.inner)
    }

    // Get a unique mutable reference to the chain if the Arc is uniquely owned. Used by `tinned_pert_multichain_insert`.
    #[inline]
    pub fn get_mut(&mut self) -> Option<&mut PertMultichain> {
        Arc::get_mut(&mut self.inner)
    }
}

// Free a perturbation multichain (NULL-safe).
#[ffi_export]
pub fn tinned_pert_multichain_free(chain: Option<PertMultichainBox>) {
    drop(chain);
}

// Clone a perturbation multichain (like Arc clone). Returns NULL on error / NULL input.
#[ffi_export]
pub fn tinned_pert_multichain_clone(
    h: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<PertMultichainBox> {
    match try_from_handle(h, "tinned_pert_multichain_clone", "PertMultichainandle", |ph| {
        PertMultichainBox::new(PertMultichainHandle::new(ph.clone_arc()))
    }) {
        Ok(b) => Some(b),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Create a new and empty perturbation multichain.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_new() -> PertMultichainBox {
    let chain = Arc::new(PertMultichain::new());
    PertMultichainBox::new(PertMultichainHandle::new(chain))
}

#[ffi_export]
pub extern "C" fn tinned_pert_multichain_from_slice(
    perturbations: PerturbationSlice<'_>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<PertMultichainBox> {
    let vec = match perturbation_vec_from_slice(perturbations, "tinned_pert_multichain_from_slice")
    {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            return None;
        },
    };

    let chain = PertMultichain::from_slice(&vec);
    Some(PertMultichainBox::new(PertMultichainHandle::new(Arc::new(chain))))
}

// Return a new PertMultichain with `p` added. Caller must free the returned box with `tinned_pert_multichain_free`.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_add(
    h: Option<&PertMultichainHandle>,
    p: Option<&PerturbationHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<PertMultichainBox> {
    // Validate + clone the perturbation Arc first.
    let pert_arc =
        match try_from_handle(p, "tinned_pert_multichain_add", "PerturbationHandle", |ph| {
            ph.clone_arc()
        }) {
            Ok(v) => v,
            Err(e) => {
                tinned_error_new(out_err, e);
                return None;
            },
        };

    // Validate the chain handle and build a new chain with the perturbation added.
    match try_with_handle(h, "tinned_pert_multichain_add", "PertMultichainHandle", |mh| {
        let new_chain = mh.as_ref().with_added_perturbation(&pert_arc);
        Ok(PertMultichainBox::new(PertMultichainHandle::new(Arc::new(new_chain))))
    }) {
        Ok(b) => Some(b),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

#[ffi_export]
pub extern "C" fn tinned_pert_multichain_insert(
    h: Option<&mut PertMultichainHandle>,
    p: Option<&PerturbationHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    // Validate + clone the perturbation Arc first.
    let pert_arc =
        match try_from_handle(p, "tinned_pert_multichain_insert", "PerturbationHandle", |ph| {
            ph.clone_arc()
        }) {
            Ok(v) => v,
            Err(e) => {
                tinned_error_new(out_err, e);
                return false;
            },
        };

    // Validate the chain handle (needs to be present and uniquely owned).
    let Some(handle) = h else {
        tinned_error_new(
            out_err,
            generic_error(
                "Null PertMultichainHandle pointer in tinned_pert_multichain_insert",
                None,
            ),
        );
        return false;
    };

    // Get a unique mutable reference to the chain (fails if shared).
    let Some(chain) = handle.get_mut() else {
        tinned_error_new(
            out_err,
            generic_error(
                "Chain handle is shared; use tinned_pert_multichain_add() or clone first",
                None,
            ),
        );
        return false;
    };

    // Perform the mutation.
    chain.insert(&pert_arc);
    true
}

// Order for a given perturbation. Returns 0 on error/`NULL` input.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_get_order(
    h: Option<&PertMultichainHandle>,
    p: Option<&PerturbationHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> u32 {
    // Validate + clone the perturbation Arc first
    let pert_arc =
        match try_from_handle(p, "tinned_pert_multichain_get_order", "PerturbationHandle", |ph| {
            ph.clone_arc()
        }) {
            Ok(v) => v,
            Err(e) => {
                tinned_error_new(out_err, e);
                return 0;
            },
        };

    match try_with_handle(h, "tinned_pert_multichain_get_order", "PertMultichainHandle", |mh| {
        Ok(mh.as_ref().get_order(&pert_arc))
    }) {
        Ok(n) => n,
        Err(e) => {
            tinned_error_new(out_err, e);
            0
        },
    }
}

// Returns `true` on error/`NULL` input.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_is_empty(
    h: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    match try_with_handle(h, "tinned_pert_multichain_is_empty", "PertMultichainHandle", |mh| {
        Ok(mh.as_ref().is_empty())
    }) {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            true
        },
    }
}

// Returns 0 on error/`NULL` input.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_total_order(
    h: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> u32 {
    match try_with_handle(h, "tinned_pert_multichain_total_order", "PertMultichainHandle", |mh| {
        Ok(mh.as_ref().total_order())
    }) {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            0
        },
    }
}

// Returns `false` on error/`NULL` input.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_is_subchain(
    lhs: Option<&PertMultichainHandle>,
    rhs: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    let res =
        try_with_handle(lhs, "tinned_pert_multichain_is_subchain", "PertMultichainHandle", |lh| {
            try_with_handle(
                rhs,
                "tinned_pert_multichain_is_subchain",
                "PertMultichainHandle",
                |rh| Ok(lh.as_ref().is_subchain(rh.as_ref())),
            )
        });

    match res {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            false
        },
    }
}

// Returns `false` on error/`NULL` input.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_is_superchain(
    lhs: Option<&PertMultichainHandle>,
    rhs: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    let res = try_with_handle(
        lhs,
        "tinned_pert_multichain_is_superchain",
        "PertMultichainHandle",
        |lh| {
            try_with_handle(
                rhs,
                "tinned_pert_multichain_is_superchain",
                "PertMultichainHandle",
                |rh| Ok(lh.as_ref().is_superchain(rh.as_ref())),
            )
        },
    );

    match res {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            false
        },
    }
}

// Returns `false` on error/`NULL` input.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_has_overlap(
    lhs: Option<&PertMultichainHandle>,
    rhs: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    let res =
        try_with_handle(lhs, "tinned_pert_multichain_has_overlap", "PertMultichainHandle", |lh| {
            try_with_handle(
                rhs,
                "tinned_pert_multichain_has_overlap",
                "PertMultichainHandle",
                |rh| Ok(lh.as_ref().has_overlap(rh.as_ref())),
            )
        });

    match res {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            false
        },
    }
}

// Display text; caller frees with `tinned_string_free`. `NULL` on error.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_display(
    h: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    match try_with_handle(h, "tinned_pert_multichain_display", "PertMultichainHandle", |mh| {
        Ok(tinned_string_to_cstr(format!("{}", mh.as_ref())))
    }) {
        Ok(s) => Some(s),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}
