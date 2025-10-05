use std::{ptr::null_mut, sync::Arc};

use tinned::expressions::AdjointMap;
use tinned::public::generic_error;

use crate::c_support::{with_downcast_expr_res, with_downcast_val};
use crate::core::{ExprBox, TinnedErrorBox, expr_box_from, set_out_err, vec_expr_from_ptrs};

#[unsafe(no_mangle)]
pub extern "C" fn tinned_adjoint_map_new(
    generator_ptrs: *const *const ExprBox,
    generator_count: usize,
    target_ptr: *const ExprBox,
    left_action: bool,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
    let generators = unsafe {
        match vec_expr_from_ptrs(generator_ptrs, generator_count, "tinned_adjoint_map_new", out_err)
        {
            Some(v) => v,
            None => return null_mut(),
        }
    };

    if target_ptr.is_null() {
        set_out_err(out_err, generic_error("Null target passed to tinned_adjoint_map_new", None));
        return null_mut();
    }

    let target = unsafe { (&*target_ptr).arc_clone() };

    match AdjointMap::new(generators, target, Some(left_action)) {
        Ok(expr) => expr_box_from(expr),
        Err(e) => {
            set_out_err(out_err, e);
            null_mut()
        },
    }
}

// Return the number of generators.
#[unsafe(no_mangle)]
pub extern "C" fn tinned_adjoint_map_generators_count(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> usize {
    with_downcast_val::<AdjointMap, usize>(h, out_err, "tinned_adjoint_map_generators_count", |a| {
        a.generators().len()
    })
}

// Return the i-th generator (cloned). Caller must unref.
#[unsafe(no_mangle)]
pub extern "C" fn tinned_adjoint_map_generator_at(
    h: *const ExprBox,
    i: usize,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
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

// Get `target` (cloned).
#[unsafe(no_mangle)]
pub extern "C" fn tinned_adjoint_map_target(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
    with_downcast_expr_res::<AdjointMap>(h, out_err, "tinned_adjoint_map_target", |a| {
        Ok(Arc::clone(a.target()))
    })
}

// Get `left_action`.
#[unsafe(no_mangle)]
pub extern "C" fn tinned_adjoint_map_left_action(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> bool {
    with_downcast_val::<AdjointMap, bool>(h, out_err, "tinned_adjoint_map_left_action", |a| {
        a.left_action()
    })
}
