use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::DotProduct;
use tinned::public::generic_error;

use crate::c_support::{ffi_map_expr_as, ffi_map_expr_as_copy};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};

#[ffi_export]
pub extern "C" fn tinned_dot_product_new(
    bra: Option<&ExprHandle>,
    use_hermitian: bool,
    ket: Option<&ExprHandle>,
    allow_braket_swap: bool,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(bra) = bra else {
        tinned_error_new(out_err, generic_error("Null bra passed to tinned_dot_product_new", None));
        return None;
    };
    let Some(ket) = ket else {
        tinned_error_new(out_err, generic_error("Null ket passed to tinned_dot_product_new", None));
        return None;
    };

    let bra_arc = bra.clone_arc();
    let ket_arc = ket.clone_arc();

    match DotProduct::new(bra_arc, use_hermitian, ket_arc, allow_braket_swap) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

impl_expr_getters!(
    DotProduct;
    tinned_dot_product_bra => |dp| Ok(Arc::clone(dp.bra())),
    tinned_dot_product_ket => |dp| Ok(Arc::clone(dp.ket())),
    tinned_dot_product_conjugate => |dp| dp.conjugate(),
);

impl_val_getters!(
    DotProduct;
    tinned_dot_product_allow_braket_swap: bool => |dp| dp.allow_braket_swap(); default = false,
);
