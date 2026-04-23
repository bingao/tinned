use safer_ffi::prelude::*;
use std::collections::BTreeMap;
use std::slice;
use std::sync::Arc;

use tinned::core::TinnedError;
use tinned::perturbations::{PertMultichain, Perturbation};
use tinned::public::generic_error;

use crate::c_support::{tinned_string_to_cstr, try_from_handle, try_with_handle};
use crate::core::{TinnedErrorBox, tinned_error_new};
use crate::perturbations::{
    PerturbationBox, PerturbationEntry, PerturbationEntrySlice, PerturbationHandle,
    PerturbationSlice, perturbation_vec_from_slice,
};

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

#[inline]
fn with_pert_multichain_ref<R>(
    h: Option<&PertMultichainHandle>,
    caller: &'static str,
    f: impl FnOnce(&PertMultichain) -> Result<R, TinnedError>,
) -> Result<R, TinnedError> {
    try_with_handle(h, caller, "PertMultichainHandle", |ph| {
        let chain = ph.as_ref();
        f(chain)
    })
}

#[inline]
fn with_two_pert_multichain_refs<R>(
    lhs: Option<&PertMultichainHandle>,
    rhs: Option<&PertMultichainHandle>,
    caller: &'static str,
    f: impl FnOnce(&PertMultichain, &PertMultichain) -> Result<R, TinnedError>,
) -> Result<R, TinnedError> {
    with_pert_multichain_ref(lhs, caller, |l| with_pert_multichain_ref(rhs, caller, |r| f(l, r)))
}

#[inline]
fn ffi_pert_multichain_return_val<R>(
    h: Option<&PertMultichainHandle>,
    caller: &'static str,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    f: impl FnOnce(&PertMultichain) -> Result<R, TinnedError>,
) -> R
where
    R: Default,
{
    match with_pert_multichain_ref(h, caller, f) {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            R::default()
        },
    }
}

#[inline]
fn ffi_pert_multichain_binary_return_val<R>(
    lhs: Option<&PertMultichainHandle>,
    rhs: Option<&PertMultichainHandle>,
    caller: &'static str,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    f: impl FnOnce(&PertMultichain, &PertMultichain) -> Result<R, TinnedError>,
) -> R
where
    R: Default,
{
    match with_two_pert_multichain_refs(lhs, rhs, caller, f) {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            R::default()
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
pub extern "C" fn tinned_pert_multichain_from_entries(
    entries: PerturbationEntrySlice,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<PertMultichainBox> {
    // Turn the raw slice into a Rust slice, with basic validation
    let raw_entries: &[PerturbationEntry] = if entries.len == 0 {
        &[]
    } else {
        if entries.ptr.is_null() {
            let err = generic_error(
                "Null PerturbationEntry slice pointer in tinned_pert_multichain_from_entries",
                None,
            );
            tinned_error_new(out_err, err);
            return None;
        }

        unsafe { slice::from_raw_parts(entries.ptr, entries.len) }
    };

    // Build the BTreeMap<Arc<Perturbation>, u32>
    let mut map: BTreeMap<Arc<Perturbation>, u32> = BTreeMap::new();

    for entry in raw_entries {
        let handle_box = entry.perturbation();
        let pert_arc = handle_box.clone_arc();
        let max_order = entry.max_order();
        // Overwrite on duplicate keys
        map.insert(pert_arc, max_order);
    }

    // Build the PertMultichain and wrap it
    let chain = PertMultichain::from_map(map);
    let handle = PertMultichainHandle::new(Arc::new(chain));
    Some(PertMultichainBox::new(handle))
}

#[ffi_export]
pub extern "C" fn tinned_pert_multichain_from_slice(
    perturbations: &PerturbationSlice,
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
    ffi_pert_multichain_return_val(h, "tinned_pert_multichain_add", out_err, |chain| {
        let new_chain = chain.with_added_perturbation(pert_arc);
        Ok(Some(PertMultichainBox::new(PertMultichainHandle::new(Arc::new(new_chain)))))
    })
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
    chain.insert(pert_arc);
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
    ffi_pert_multichain_return_val(h, "tinned_pert_multichain_get_order", out_err, |chain| {
        Ok(chain.get_order(&pert_arc))
    })
}

// Returns `true` on error/`NULL` input.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_is_empty(
    h: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    ffi_pert_multichain_return_val(h, "tinned_pert_multichain_is_empty", out_err, |chain| {
        Ok(chain.is_empty())
    })
}

// Returns 0 on error/`NULL` input.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_total_order(
    h: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> u32 {
    ffi_pert_multichain_return_val(h, "tinned_pert_multichain_total_order", out_err, |chain| {
        Ok(chain.total_order())
    })
}

// Returns a cloned vector of the perturbation multichain
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_to_vec(
    h: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> repr_c::Vec<PerturbationBox> {
    match try_with_handle(h, "tinned_pert_multichain_to_vec", "PertMultichainHandle", |ph| {
        let vec_perts = ph.as_ref().to_vec();

        // Build a standard Vec first
        let mut vec_pert_boxes: Vec<PerturbationBox> = Vec::with_capacity(vec_perts.len());
        for pert in vec_perts {
            vec_pert_boxes.push(PerturbationBox::new(PerturbationHandle::new(pert)));
        }

        // Then convert to repr_c::Vec
        let out: repr_c::Vec<PerturbationBox> = vec_pert_boxes.into();
        Ok(out)
    }) {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            // Build an empty std Vec and convert to repr_c::Vec as the fallback
            Vec::<PerturbationBox>::new().into()
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
    ffi_pert_multichain_binary_return_val(
        lhs,
        rhs,
        "tinned_pert_multichain_is_subchain",
        out_err,
        |l, r| Ok(l.is_subchain(r)),
    )
}

// Returns `false` on error/`NULL` input.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_is_superchain(
    lhs: Option<&PertMultichainHandle>,
    rhs: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    ffi_pert_multichain_binary_return_val(
        lhs,
        rhs,
        "tinned_pert_multichain_is_superchain",
        out_err,
        |l, r| Ok(l.is_superchain(r)),
    )
}

// Returns `false` on error/`NULL` input.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_has_overlap(
    lhs: Option<&PertMultichainHandle>,
    rhs: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    ffi_pert_multichain_binary_return_val(
        lhs,
        rhs,
        "tinned_pert_multichain_has_overlap",
        out_err,
        |l, r| Ok(l.has_overlap(r)),
    )
}

// Display text; caller frees with `tinned_string_free`. `NULL` on error.
#[ffi_export]
pub extern "C" fn tinned_pert_multichain_display(
    h: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    ffi_pert_multichain_return_val(h, "tinned_pert_multichain_display", out_err, |chain| {
        Ok(Some(tinned_string_to_cstr(format!("{}", chain))))
    })
}
