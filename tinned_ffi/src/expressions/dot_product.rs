use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::DotProduct;
use tinned::public::generic_error;

use crate::c_support::{with_downcast_expr_res, with_downcast_val};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, set_out_err};

#[ffi_export]
pub extern "C" fn tinned_dot_product_new(
    bra: Option<&ExprHandle>,
    use_hermitian: bool,
    ket: Option<&ExprHandle>,
    allow_braket_swap: bool,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(bra) = bra else {
        set_out_err(out_err, generic_error("Null bra passed to tinned_dot_product_new", None));
        return None;
    };
    let Some(ket) = ket else {
        set_out_err(out_err, generic_error("Null ket passed to tinned_dot_product_new", None));
        return None;
    };

    let bra_arc = bra.clone_arc();
    let ket_arc = ket.clone_arc();

    match DotProduct::new(bra_arc, use_hermitian, ket_arc, allow_braket_swap) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            set_out_err(out_err, e);
            None
        },
    }
}

// Get `bra` (cloned).
#[ffi_export]
pub extern "C" fn tinned_dot_product_bra(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    with_downcast_expr_res::<DotProduct>(h, out_err, "tinned_dot_product_bra", |dp| {
        Ok(Arc::clone(dp.bra()))
    })
}

// Get `ket` (cloned).
#[ffi_export]
pub extern "C" fn tinned_dot_product_ket(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    with_downcast_expr_res::<DotProduct>(h, out_err, "tinned_dot_product_ket", |dp| {
        Ok(Arc::clone(dp.ket()))
    })
}

// Get `allow_braket_swap`.
#[ffi_export]
pub extern "C" fn tinned_dot_product_allow_braket_swap(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    with_downcast_val::<DotProduct, bool>(
        h,
        out_err,
        "tinned_dot_product_allow_braket_swap",
        |dp| dp.allow_braket_swap(),
    )
    .unwrap_or(false)
}

// Compute `conjugate()` and return a new expression.
#[ffi_export]
pub extern "C" fn tinned_dot_product_conjugate(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    with_downcast_expr_res::<DotProduct>(h, out_err, "tinned_dot_product_conjugate", |dp| {
        dp.conjugate()
    })
}
