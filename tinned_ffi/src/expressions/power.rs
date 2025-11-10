use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::Power;
use tinned::public::generic_error;

use crate::c_support::{with_downcast_expr, with_downcast_val};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};

#[ffi_export]
pub extern "C" fn tinned_power_new(
    base: Option<&ExprHandle>,
    exponent: i64,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(base) = base else {
        tinned_error_new(out_err, generic_error("Null base passed to tinned_power_new", None));
        return None;
    };
    let base_arc = base.clone_arc();

    match Power::new(base_arc, exponent) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Get `base` (cloned).
impl_expr_getters!(
    Power;
    tinned_power_base => |p| Ok(Arc::clone(p.base())),
);

// Get `exponent`.
impl_val_getters!(
    Power;
    tinned_power_exponent: i64 => |p| p.exponent(); default = 0,
);
