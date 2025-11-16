use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::HermitianTranspose;
use tinned::public::generic_error;

use crate::c_support::ffi_map_expr_as;
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};

#[ffi_export]
pub extern "C" fn tinned_hermitian_transpose_new(
    argument: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(argument) = argument else {
        tinned_error_new(
            out_err,
            generic_error("Null argument passed to tinned_hermitian_transpose_new", None),
        );
        return None;
    };
    let arg_arc = argument.clone_arc();

    match HermitianTranspose::new(arg_arc) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Get `argument` (cloned).
impl_expr_getters!(
    HermitianTranspose;
    tinned_hermitian_transpose_argument => |herm| Ok(Arc::clone(herm.argument())),
);
