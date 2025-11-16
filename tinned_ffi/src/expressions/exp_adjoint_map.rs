use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::ExpAdjointMap;
use tinned::public::generic_error;

use crate::c_support::{ffi_map_expr_as, ffi_map_expr_as_copy};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::{PertMultichainHandle, PertMultichainBox};

#[ffi_export]
pub extern "C" fn tinned_exp_adjoint_map_new(
    generator: Option<&ExprHandle>,
    target: Option<&ExprHandle>,
    left_action: bool,
    max_fold: u32,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(generator) = generator else {
        tinned_error_new(
            out_err,
            generic_error("Null generator passed to tinned_exp_adjoint_map_new", None),
        );
        return None;
    };
    let generator_arc = generator.clone_arc();

    let Some(target) = target else {
        tinned_error_new(
            out_err,
            generic_error("Null target passed to tinned_exp_adjoint_map_new", None),
        );
        return None;
    };
    let target_arc = target.clone_arc();

    match ExpAdjointMap::builder(generator_arc, target_arc)
        .left_action(left_action)
        .max_fold(max_fold)
        .build()
    {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

impl_expr_getters!(
    ExpAdjointMap;
    tinned_exp_adjoint_map_generator => |ead| Ok(Arc::clone(ead.generator())),
    tinned_exp_adjoint_map_target => |ead| Ok(Arc::clone(ead.target())),
    tinned_exp_adjoint_map_result => |ead| Ok(Arc::clone(ead.result())),
);

impl_val_getters!(
    ExpAdjointMap;
    tinned_exp_adjoint_map_is_temporum: bool => |ead| ead.is_temporum(); default = false,
    tinned_exp_adjoint_map_left_action: bool => |ead| ead.left_action(); default = false,
    tinned_exp_adjoint_map_max_fold: u32 => |ead| ead.max_fold(); default = 0,
    tinned_exp_adjoint_map_is_zero_strength: bool => |ead| ead.is_zero_strength(); default = false,
);

// Get `derivative` (cloned).
impl_pert_multichain_getter!(
    tinned_exp_adjoint_map_derivative : ExpAdjointMap => |ead| ead.derivative().clone()
);
