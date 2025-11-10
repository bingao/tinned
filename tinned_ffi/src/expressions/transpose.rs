use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::Transpose;
use tinned::public::generic_error;

use crate::c_support::with_downcast_expr;
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};

#[ffi_export]
pub extern "C" fn tinned_transpose_new(
    argument: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(argument) = argument else {
        tinned_error_new(
            out_err,
            generic_error("Null argument passed to tinned_transpose_new", None),
        );
        return None;
    };
    let arg_arc = argument.clone_arc();

    match Transpose::new(arg_arc) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Get `argument` (cloned).
impl_expr_getters!(
    Transpose;
    tinned_transpose_argument => |tr| Ok(Arc::clone(tr.argument())),
);
