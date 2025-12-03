use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::core::Expr;
use tinned::expressions::Symbol;
use tinned::public::generic_error;

use crate::c_support::{ffi_map_expr_as, tinned_string_from_cstr, tinned_string_to_cstr};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};

/// Create a new `Symbol` expression.
/// - `name`: UTF-8 C string (nullable). On NULL/invalid, sets `out_err` and returns `None`.
/// - Returns an `ExprBox` on success.
#[ffi_export]
pub extern "C" fn tinned_symbol_new(
    name: Option<char_p::Ref<'_>>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(name) = tinned_string_from_cstr(name) else {
        tinned_error_new(
            out_err,
            generic_error("Null or invalid name passed to tinned_symbol_new", None),
        );
        return None;
    };

    // `Symbol::new` is infallible and returns Arc<dyn Expr>.
    let expr_arc: Arc<dyn Expr> = Symbol::new(name);
    Some(ExprBox::new(ExprHandle::new(expr_arc)))
}

// Get the `name` of a `Symbol`.
// - Returns a newly allocated C string; free with `tinned_string_free`.
impl_cstr_getter!(
    tinned_symbol_name : Symbol => |s| s.name().to_string()
);
