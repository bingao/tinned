use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::Conjugate;
use tinned::public::generic_error;

use crate::c_support::with_downcast_expr_res;
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, set_out_err};

#[ffi_export]
pub extern "C" fn tinned_conjugate_new(
    argument: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(argument) = argument else {
        set_out_err(out_err, generic_error("Null argument passed to tinned_conjugate_new", None));
        return None;
    };
    let arg_arc = argument.clone_arc();

    match Conjugate::new(arg_arc) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            set_out_err(out_err, e);
            None
        },
    }
}

// Get `argument` (cloned).
#[ffi_export]
pub extern "C" fn tinned_conjugate_argument(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    with_downcast_expr_res::<Conjugate>(h, out_err, "tinned_conjugate_argument", |cj| {
        Ok(Arc::clone(cj.argument()))
    })
}
