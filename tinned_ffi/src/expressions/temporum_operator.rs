use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::TemporumOperator;
use tinned::public::generic_error;

use crate::c_support::{with_downcast_expr, with_downcast_pert_multichain, with_downcast_val};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::PertMultichainBox;

#[ffi_export]
pub extern "C" fn tinned_temporum_operator_new(
    argument: Option<&ExprHandle>,
    is_forward: bool,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(argument) = argument else {
        tinned_error_new(
            out_err,
            generic_error("Null argument passed to tinned_temporum_operator_new", None),
        );
        return None;
    };
    let argument_arc = argument.clone_arc();

    match TemporumOperator::builder(argument_arc).is_forward(is_forward).build() {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

impl_val_getters!(
    TemporumOperator;
    tinned_temporum_operator_is_forward: bool => |op| op.is_forward(); default = false,
);

impl_expr_getters!(
    TemporumOperator;
    tinned_temporum_operator_argument => |op| Ok(Arc::clone(op.argument())),
    tinned_temporum_operator_frequency => |op| op.frequency(),
);

// Get `derivative` (cloned).
impl_pert_multichain_getter!(
    tinned_temporum_operator_derivative : TemporumOperator => |op| op.derivative().map(|mc| Arc::new(mc.clone()))
);
