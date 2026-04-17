use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::TimeEvolution;
use tinned::public::generic_error;

use crate::c_support::ffi_map_expr_as;
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::{PertMultichainBox, PertMultichainHandle};

#[ffi_export]
pub extern "C" fn tinned_time_evolution_new(
    argument: Option<&ExprHandle>,
    is_forward: bool,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(argument) = argument else {
        tinned_error_new(
            out_err,
            generic_error("Null argument passed to tinned_time_evolution_new", None),
        );
        return None;
    };
    let argument_arc = argument.clone_arc();

    match TimeEvolution::builder(argument_arc).is_forward(is_forward).build() {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

impl_val_getters!(
    TimeEvolution;
    tinned_time_evolution_is_forward: bool => |op| op.is_forward(); default = false,
);

impl_expr_getters!(
    TimeEvolution;
    tinned_time_evolution_argument => |op| Ok(Arc::clone(op.argument())),
    tinned_time_evolution_frequency => |op| op.frequency(),
);

// Get `derivative` (cloned).
#[ffi_export]
pub extern "C" fn tinned_time_evolution_derivative(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<PertMultichainBox> {
    ffi_map_expr_as::<TimeEvolution, _>(h, out_err, "tinned_time_evolution_derivative", |op| {
        op.derivative()
            .map(|chain| PertMultichainBox::new(PertMultichainHandle::new(Arc::new(chain.clone()))))
    })
}
