use safer_ffi::prelude::*;

use tinned::public::NumberTolerance;

#[derive_ReprC]
#[repr(opaque)]
pub struct NumberToleranceHandle {
    inner: NumberTolerance,
}

pub type NumberToleranceBox = repr_c::Box<NumberToleranceHandle>;

impl NumberToleranceHandle {
    #[inline]
    pub fn new(tol: NumberTolerance) -> Self {
        Self {
            inner: tol,
        }
    }

    #[inline]
    pub fn as_ref(&self) -> &NumberTolerance {
        &self.inner
    }

    //#[inline]
    //pub fn into_inner(self: repr_c::Box<Self>) -> NumberTolerance { repr_c::Box::into_inner(self).inner }
}

#[ffi_export]
pub fn tinned_number_tolerance_new(abs_error: f64, rel_error: f64) -> Option<NumberToleranceBox> {
    // If the constructor panics on invalid values, prevent unwinding across FFI:
    let tol = std::panic::catch_unwind(move || NumberTolerance::new(abs_error, rel_error));
    match tol {
        Ok(val) => Some(NumberToleranceBox::new(NumberToleranceHandle::new(val))),
        Err(_) => None,
    }
}

#[ffi_export]
pub fn tinned_number_tolerance_zero() -> NumberToleranceBox {
    NumberToleranceBox::new(NumberToleranceHandle::new(NumberTolerance::zero()))
}

#[ffi_export]
pub fn tinned_number_tolerance_abs_error(h: &NumberToleranceHandle) -> f64 {
    h.as_ref().abs_error()
}

#[ffi_export]
pub fn tinned_number_tolerance_rel_error(h: &NumberToleranceHandle) -> f64 {
    h.as_ref().rel_error()
}

#[ffi_export]
pub fn tinned_number_tolerance_free(tol: Option<NumberToleranceBox>) {
    drop(tol)
}

#[ffi_export]
pub fn tinned_get_global_number_tolerance() -> NumberToleranceBox {
    let copy = tinned::get_number_tolerance();
    NumberToleranceBox::new(NumberToleranceHandle::new(copy))
}

#[ffi_export]
pub fn tinned_set_global_number_tolerance(h: &NumberToleranceHandle) {
    tinned::set_number_tolerance(h.as_ref().clone());
}
