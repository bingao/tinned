use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::AdjointMap;
use tinned::public::generic_error;

use crate::c_support::{with_downcast_expr_res, with_downcast_val};
use crate::core::{
    ExprBox, ExprHandle, ExprSlice, TinnedErrorBox, expr_vec_from_slice, set_out_err,
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
            set_out_err(out_err, e);
            return None;
        },
    };

    let Some(target) = target else {
        set_out_err(out_err, generic_error("Null target passed to tinned_adjoint_map_new", None));
        return None;
    };
    let target_arc = target.clone_arc();

    match AdjointMap::new(generators_vec, target_arc, Some(left_action)) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            set_out_err(out_err, e);
            None
        },
    }
}

#[ffi_export]
pub extern "C" fn tinned_adjoint_map_generators_count(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> usize {
    with_downcast_val::<AdjointMap, usize>(h, out_err, "tinned_adjoint_map_generators_count", |a| {
        a.generators().len()
    })
    .unwrap_or(0)
}

// Return the i-th generator (cloned). Caller must free the returned ExprBox.
#[ffi_export]
pub extern "C" fn tinned_adjoint_map_generator_at(
    h: Option<&ExprHandle>,
    i: usize,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    with_downcast_expr_res::<AdjointMap>(h, out_err, "tinned_adjoint_map_generator_at", |a| {
        a.generators().get(i).cloned().ok_or_else(|| {
            generic_error(
                format!(
                    "Index {i} out of bounds (len = {}) in tinned_adjoint_map_generator_at",
                    a.generators().len()
                ),
                None,
            )
        })
    })
}

#[ffi_export]
pub extern "C" fn tinned_adjoint_map_target(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    with_downcast_expr_res::<AdjointMap>(h, out_err, "tinned_adjoint_map_target", |a| {
        Ok(Arc::clone(a.target()))
    })
}

#[ffi_export]
pub extern "C" fn tinned_adjoint_map_left_action(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    with_downcast_val::<AdjointMap, bool>(h, out_err, "tinned_adjoint_map_left_action", |a| {
        a.left_action()
    })
    .unwrap_or(false)
}
