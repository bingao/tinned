use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::AdjointMap;
use tinned::public::generic_error;

use crate::c_support::{with_downcast_expr, with_downcast_val};
use crate::core::{
    ExprBox, ExprHandle, ExprSlice, TinnedErrorBox, expr_vec_from_slice, tinned_error_new,
};

#[ffi_export]
pub extern "C" fn tinned_adjoint_map_new(
    generators: ExprSlice<'_>,
    target: Option<&ExprHandle>,
    left_action: bool,
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

    match AdjointMap::new(generators_vec, target_arc, Some(left_action)) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

impl_val_getters!(
    AdjointMap;
    tinned_adjoint_map_generators_count: usize => |a| a.generators().len(); default = 0,
    tinned_adjoint_map_left_action: bool => |a| a.left_action(); default = false,
);

// Return the i-th generator (cloned). Caller must free the returned ExprBox.
impl_expr_index_getter!(tinned_adjoint_map_generator_at : AdjointMap => generators);

impl_expr_getters!(
    AdjointMap;
    tinned_adjoint_map_target => |adj| Ok(Arc::clone(adj.target())),
);
