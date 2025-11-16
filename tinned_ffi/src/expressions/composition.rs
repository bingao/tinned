use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::Composition;
use tinned::public::generic_error;

use crate::c_support::{
    tinned_string_from_cstr, ffi_map_expr_as_copy, ffi_map_expr_as, tinned_string_to_cstr,
};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};

#[ffi_export]
pub extern "C" fn tinned_composition_new(
    name: Option<char_p::Ref<'_>>,
    order: u32,
    inner: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(name) = tinned_string_from_cstr(name) else {
        tinned_error_new(
            out_err,
            generic_error("Null or invalid name passed to tinned_composition_new", None),
        );
        return None;
    };

    let Some(inner) = inner else {
        tinned_error_new(
            out_err,
            generic_error("Null inner passed to tinned_composition_new", None),
        );
        return None;
    };
    let inner_arc = inner.clone_arc();

    match Composition::new(name, order, inner_arc) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Get `name` (caller must free the returned C string).
impl_cstr_getter!(
    tinned_composition_name : Composition => |c| c.name().to_string()
);

// Get `order`.
impl_val_getters!(
    Composition;
    tinned_composition_order: u32 => |c| c.order(); default = 0,
);

// Get `inner` (cloned).
impl_expr_getters!(
    Composition;
    tinned_composition_inner => |c| Ok(Arc::clone(c.inner())),
);
