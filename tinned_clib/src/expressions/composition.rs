use std::{os::raw::c_char, ptr::null_mut, sync::Arc};

use tinned::expressions::Composition;
use tinned::public::generic_error;

use crate::c_support::{
    cstr_to_string, with_downcast_cstr, with_downcast_expr_res, with_downcast_val,
};
use crate::core::{ExprBox, TinnedErrorBox, expr_box_from, set_out_err};

#[unsafe(no_mangle)]
pub extern "C" fn tinned_composition_new(
    name_cstr: *const c_char,
    order: u32,
    inner_ptr: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
    let Some(name) = cstr_to_string(name_cstr) else {
        set_out_err(
            out_err,
            generic_error("Null or invalid name passed to tinned_composition_new", None),
        );
        return null_mut();
    };
    if inner_ptr.is_null() {
        set_out_err(out_err, generic_error("Null inner passed to tinned_composition_new", None));
        return null_mut();
    }
    let inner = unsafe { (&*inner_ptr).arc_clone() };

    match Composition::new(name, order, inner) {
        Ok(expr) => expr_box_from(expr),
        Err(e) => {
            set_out_err(out_err, e);
            null_mut()
        },
    }
}

// Get `name` (C must free).
#[unsafe(no_mangle)]
pub extern "C" fn tinned_composition_name(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut c_char {
    with_downcast_cstr::<Composition>(h, out_err, "tinned_composition_name", |c| {
        c.name().to_string()
    })
}

// Get `order`.
#[unsafe(no_mangle)]
pub extern "C" fn tinned_composition_order(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> u32 {
    with_downcast_val::<Composition, u32>(h, out_err, "tinned_composition_order", |c| c.order())
}

// Get `inner` (cloned).
#[unsafe(no_mangle)]
pub extern "C" fn tinned_composition_inner(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
    with_downcast_expr_res::<Composition>(h, out_err, "tinned_composition_inner", |c| {
        Ok(Arc::clone(c.inner()))
    })
}
