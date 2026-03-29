use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::{AdjointMap, AdjointMode};
use tinned::public::generic_error;

use crate::c_support::ffi_map_expr_as_exprvec;
use crate::core::{
    ExprBox, ExprHandle, ExprSlice, TinnedErrorBox, expr_vec_from_slice, tinned_error_new,
};

#[ffi_export]
pub extern "C" fn tinned_adjoint_map_new(
    generators: &ExprSlice,
    target: Option<&ExprHandle>,
    left_action: bool,
    adjoint_mode: AdjointMode,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let generators_vec = match expr_vec_from_slice(generators, "tinned_adjoint_map_new") {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            return None;
        },
    };

    let Some(target) = target else {
        tinned_error_new(
            out_err,
            generic_error("Null target passed to tinned_adjoint_map_new", None),
        );
        return None;
    };
    let target_arc = target.clone_arc();

    match AdjointMap::new(generators_vec, target_arc, Some(left_action), Some(adjoint_mode)) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

impl_val_getters!(
    AdjointMap;
    tinned_adjoint_map_left_action: bool => |a| a.left_action(); default = false,
    tinned_adjoint_map_adjoint_mode: AdjointMode => |a| a.adjoint_mode(); default = AdjointMode::Commutative,
);

// Returns a cloned vector of generators
#[ffi_export]
pub extern "C" fn tinned_adjoint_map_generators(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> repr_c::Vec<ExprBox> {
    ffi_map_expr_as_exprvec::<AdjointMap>(h, out_err, "tinned_adjoint_map_generators", |adj| {
        adj.generators()
    })
}

impl_expr_getters!(
    AdjointMap;
    tinned_adjoint_map_target => |adj| Ok(Arc::clone(adj.target())),
);
