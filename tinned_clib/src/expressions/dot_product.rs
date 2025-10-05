use std::{ptr::null_mut, sync::Arc};

use tinned::expressions::DotProduct;
use tinned::public::generic_error;

use crate::c_support::{with_downcast_expr_res, with_downcast_val};
use crate::core::{ExprBox, TinnedErrorBox, expr_box_from, set_out_err};

#[unsafe(no_mangle)]
pub extern "C" fn tinned_dot_product_new(
    bra_ptr: *const ExprBox,
    use_hermitian: bool,
    ket_ptr: *const ExprBox,
    allow_braket_swap: bool,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
    if bra_ptr.is_null() || ket_ptr.is_null() {
        set_out_err(out_err, generic_error("Null bra/ket passed to tinned_dot_product_new", None));
        return null_mut();
    }
    let bra = unsafe { (&*bra_ptr).arc_clone() };
    let ket = unsafe { (&*ket_ptr).arc_clone() };

    match DotProduct::new(bra, use_hermitian, ket, allow_braket_swap) {
        Ok(expr) => expr_box_from(expr),
        Err(e) => {
            set_out_err(out_err, e);
            null_mut()
        },
    }
}

// Get `bra` (cloned).
#[unsafe(no_mangle)]
pub extern "C" fn tinned_dot_product_bra(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
    with_downcast_expr_res::<DotProduct>(h, out_err, "tinned_dot_product_bra", |dp| {
        Ok(Arc::clone(dp.bra()))
    })
}

// Get `ket` (cloned).
#[unsafe(no_mangle)]
pub extern "C" fn tinned_dot_product_ket(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
    with_downcast_expr_res::<DotProduct>(h, out_err, "tinned_dot_product_ket", |dp| {
        Ok(Arc::clone(dp.ket()))
    })
}

// Get `allow_braket_swap`.
#[unsafe(no_mangle)]
pub extern "C" fn tinned_dot_product_allow_braket_swap(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> bool {
    with_downcast_val::<DotProduct, bool>(
        h,
        out_err,
        "tinned_dot_product_allow_braket_swap",
        |dp| dp.allow_braket_swap(),
    )
}

// Compute `conjugate()` and return a new expression.
#[unsafe(no_mangle)]
pub extern "C" fn tinned_dot_product_conjugate(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
    with_downcast_expr_res::<DotProduct>(h, out_err, "tinned_dot_product_conjugate", |dp| {
        dp.conjugate()
    })
}
