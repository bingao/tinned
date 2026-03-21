use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::core::Expr;
use tinned::expressions::ExcitationOperator;
use tinned::public::generic_error;

use crate::c_support::tinned_string_from_cstr;
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};

/// Create a new `ExcitationOperator` expression.
/// - `name`: UTF-8 C string (nullable). On NULL/invalid, sets `out_err` and returns `None`.
/// - Returns an `ExprBox` on success.
#[ffi_export]
pub extern "C" fn tinned_excitation_operator_new(
    name: Option<char_p::Ref<'_>>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(name) = tinned_string_from_cstr(name) else {
        tinned_error_new(
            out_err,
            generic_error("Null or invalid name passed to tinned_excitation_operator_new", None),
        );
        return None;
    };

    // `ExcitationOperator::new` is infallible and returns Arc<dyn Expr>.
    let expr_arc: Arc<dyn Expr> = ExcitationOperator::new(name);
    Some(ExprBox::new(ExprHandle::new(expr_arc)))
}

// Get the `name` of an `ExcitationOperator`.
// - Returns a newly allocated C string; free with `tinned_string_free`.
impl_cstr_getter!(
    tinned_excitation_operator_name: ExcitationOperator => |s| s.name().to_string()
);
