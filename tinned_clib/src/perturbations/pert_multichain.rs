use std::{os::raw::c_char, ptr::null_mut, sync::Arc};

use tinned::perturbations::{PertMultichain, Perturbation};
use tinned::public::generic_error;

use crate::c_support::to_cstring;
use crate::core::{TinnedErrorBox, set_out_err};
use crate::perturbations::{
    PertMultichainBox, PerturbationBox, vec_pert_from_ptrs, with_pert_multichain_or_err,
};

#[unsafe(no_mangle)]
pub extern "C" fn tinned_pert_multichain_ref(h: *mut PertMultichainBox) -> *mut PertMultichainBox {
    if h.is_null() {
        return null_mut();
    }
    let chain = unsafe { &*h }.arc_clone();
    PertMultichainBox::new(chain).into_raw()
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_pert_multichain_unref(h: *mut PertMultichainBox) {
    if h.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(h));
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_pert_multichain_new() -> *mut PertMultichainBox {
    PertMultichainBox::new(PertMultichain::new()).into_raw()
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_pert_multichain_from_slice(
    perturbation_ptrs: *const *const PerturbationBox,
    perturbation_count: usize,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut PertMultichainBox {
    let slice: Vec<Arc<Perturbation>> = unsafe {
        match vec_pert_from_ptrs(
            perturbation_ptrs,
            perturbation_count,
            "tinned_pert_multichain_from_slice",
            out_err,
        ) {
            Some(v) => v,
            None => return null_mut(),
        }
    };

    let chain = PertMultichain::from_slice(&slice);
    PertMultichainBox::new(chain).into_raw()
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_pert_multichain_insert(
    h: *mut PertMultichainBox,
    p: *const PerturbationBox,
    out_err: *mut *mut TinnedErrorBox,
) -> bool {
    if h.is_null() || p.is_null() {
        set_out_err(
            out_err,
            generic_error(
                "Null pointer(s) to chain or perturbation passed to tinned_pert_multichain_insert",
                None,
            ),
        );
        return false;
    }
    let pert = unsafe { &*p }.arc_clone();

    // We store `Arc<PertMultichain>` inside the box. To call `&mut self` methods,
    // we must have unique ownership of that `Arc`.
    let chain_mut: &mut Arc<PertMultichain> = unsafe { &mut *h }.as_arc_mut();
    match Arc::get_mut(chain_mut) {
        Some(chain) => {
            chain.insert(&pert);
            true
        },
        None => {
            // Not uniquely owned. Suggest using the functional add instead.
            set_out_err(
                out_err,
                generic_error(
                    "Chain handle is shared; use tinned_pert_multichain_add() or clone first",
                    None,
                ),
            );
            false
        },
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_pert_multichain_get_order(
    h: *const PertMultichainBox,
    p: *const PerturbationBox,
) -> u32 {
    if h.is_null() || p.is_null() {
        return 0;
    }
    let pert = unsafe { &*p }.arc_clone();
    unsafe { &*h }.as_ref().get_order(&pert)
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_pert_multichain_is_empty(h: *const PertMultichainBox) -> bool {
    if h.is_null() {
        return true;
    }
    unsafe { &*h }.as_ref().is_empty()
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_pert_multichain_total_order(h: *const PertMultichainBox) -> u32 {
    if h.is_null() {
        return 0;
    }
    unsafe { &*h }.as_ref().total_order()
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_pert_multichain_is_subchain(
    lhs: *const PertMultichainBox,
    rhs: *const PertMultichainBox,
) -> bool {
    if lhs.is_null() || rhs.is_null() {
        return false;
    }
    unsafe { &*lhs }.as_ref().is_subchain(unsafe { &*rhs }.as_ref())
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_pert_multichain_is_superchain(
    lhs: *const PertMultichainBox,
    rhs: *const PertMultichainBox,
) -> bool {
    if lhs.is_null() || rhs.is_null() {
        return false;
    }
    unsafe { &*lhs }.as_ref().is_superchain(unsafe { &*rhs }.as_ref())
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_pert_multichain_has_overlap(
    lhs: *const PertMultichainBox,
    rhs: *const PertMultichainBox,
) -> bool {
    if lhs.is_null() || rhs.is_null() {
        return false;
    }
    unsafe { &*lhs }.as_ref().has_overlap(unsafe { &*rhs }.as_ref())
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_pert_multichain_display(
    h: *const PertMultichainBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut c_char {
    with_pert_multichain_or_err(h, out_err, "tinned_pert_multichain_display", |chain| {
        format!("{}", chain)
    })
    .map(to_cstring)
    .unwrap_or(null_mut())
}
