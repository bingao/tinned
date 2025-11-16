use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::TemporumOverlap;
use tinned::perturbations::PertMultichain;
use tinned::public::generic_error;

use crate::c_support::{ffi_map_expr_as, ffi_map_expr_as_copy};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::{PertMultichainBox, PertMultichainHandle};

#[ffi_export]
pub extern "C" fn tinned_temporum_overlap_new(
    dependencies: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(dependencies) = dependencies else {
        tinned_error_new(
            out_err,
            generic_error("Null dependencies passed to tinned_temporum_overlap_new", None),
        );
        return None;
    };
    let deps: PertMultichain = dependencies.as_ref().clone();

    match TemporumOverlap::builder(deps).build() {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

impl_val_getters!(
    TemporumOverlap;
    tinned_temporum_overlap_is_zero_strength: bool => |op| op.is_zero_strength(); default = false,
);

impl_expr_getters!(
    TemporumOverlap;
    tinned_temporum_overlap_braket => |op| Ok(Arc::clone(op.braket())),
);

// Get `derivative` (cloned).
impl_pert_multichain_getter!(
    tinned_temporum_overlap_derivative : TemporumOverlap => |op| op.derivative().clone()
);

// Get `dependencies` (cloned).
impl_pert_multichain_getter!(
    tinned_temporum_overlap_dependencies : TemporumOverlap => |op| op.dependencies().clone()
);
