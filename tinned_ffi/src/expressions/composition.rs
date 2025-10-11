use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::Composition;
use tinned::public::generic_error;

use crate::c_support::{
    tinned_string_from_cstr, with_downcast_cstr, with_downcast_expr_res, with_downcast_val,
};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, set_out_err};

#[ffi_export]
pub extern "C" fn tinned_composition_new(
    name: Option<char_p::Ref<'_>>,
    order: u32,
    inner: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(name) = tinned_string_from_cstr(name) else {
        set_out_err(
            out_err,
            generic_error("Null or invalid name passed to tinned_composition_new", None),
        );
        return None;
    };

    let Some(inner) = inner else {
        set_out_err(out_err, generic_error("Null inner passed to tinned_composition_new", None));
        return None;
    };
    let inner_arc = inner.clone_arc();

    match Composition::new(name, order, inner_arc) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            set_out_err(out_err, e);
            None
        },
    }
}

// Get `name` (caller must free the returned C string).
#[ffi_export]
pub extern "C" fn tinned_composition_name(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    with_downcast_cstr::<Composition>(h, out_err, "tinned_composition_name", |c| {
        c.name().to_string()
    })
}

// Get `order`.
#[ffi_export]
pub extern "C" fn tinned_composition_order(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> u32 {
    with_downcast_val::<Composition, u32>(h, out_err, "tinned_composition_order", |c| c.order())
        .unwrap_or(0)
}

// Get `inner` (cloned).
#[ffi_export]
pub extern "C" fn tinned_composition_inner(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    with_downcast_expr_res::<Composition>(h, out_err, "tinned_composition_inner", |c| {
        Ok(Arc::clone(c.inner()))
    })
}
