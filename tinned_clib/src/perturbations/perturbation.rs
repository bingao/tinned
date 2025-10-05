use std::{os::raw::c_char, ptr::null_mut, sync::Arc};

use tinned::perturbations::Perturbation;
use tinned::public::generic_error;

use crate::c_support::{cstr_to_string, to_cstring};
use crate::core::{ExprBox, TinnedErrorBox, expr_box_from, set_out_err};
use crate::perturbations::{PerturbationBox, perturbation_box_from, with_perturbation_or_err};

#[unsafe(no_mangle)]
pub extern "C" fn tinned_perturbation_ref(h: *mut PerturbationBox) -> *mut PerturbationBox {
    if h.is_null() {
        return null_mut();
    }
    let pert = unsafe { &*h }.arc_clone();
    perturbation_box_from(pert)
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_perturbation_unref(h: *mut PerturbationBox) {
    if h.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(h));
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_perturbation_new(
    name_cstr: *const c_char,
    frequency_ptr: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut PerturbationBox {
    let Some(name) = cstr_to_string(name_cstr) else {
        set_out_err(
            out_err,
            generic_error("Null perturbation name passed to tinned_perturbation_new", None),
        );
        return null_mut();
    };
    if frequency_ptr.is_null() {
        set_out_err(
            out_err,
            generic_error("Null frequency passed to tinned_perturbation_new", None),
        );
        return null_mut();
    }
    let frequency = unsafe { &*frequency_ptr }.arc_clone();

    let pert = Perturbation::new(name, frequency);
    perturbation_box_from(pert)
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_perturbation_name(
    h: *const PerturbationBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut c_char {
    with_perturbation_or_err(h, out_err, "tinned_perturbation_name", |p| to_cstring(p.name()))
        .unwrap_or(null_mut())
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_perturbation_frequency(
    h: *const PerturbationBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
    with_perturbation_or_err(h, out_err, "tinned_perturbation_frequency", |p| {
        Arc::clone(p.frequency())
    })
    .map(expr_box_from)
    .unwrap_or(null_mut())
}

#[unsafe(no_mangle)]
pub extern "C" fn tinned_perturbation_display(
    h: *const PerturbationBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut c_char {
    with_perturbation_or_err(h, out_err, "tinned_perturbation_display", |p| format!("{}", p))
        .map(to_cstring)
        .unwrap_or(null_mut())
}
