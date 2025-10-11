use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::core::Expr;
use tinned::public::generic_error;

use crate::c_support::{
    tinned_string_from_cstr, tinned_string_to_cstr, try_from_handle, try_with_handle,
};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, set_out_err};

// Free an expression (NULL-safe).
#[ffi_export]
pub fn tinned_expr_free(expr: Option<ExprBox>) {
    drop(expr);
}

// Clone an expression (like Arc clone). Returns NULL on error / NULL input.
#[ffi_export]
pub fn tinned_expr_clone(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    match try_from_handle(h, "tinned_expr_clone", "ExprHandle", |eh| {
        ExprBox::new(ExprHandle::new(eh.clone_arc()))
    }) {
        Ok(b) => Some(b),
        Err(e) => {
            set_out_err(out_err, e);
            None
        },
    }
}

// Whether the expression is scalar. Returns false on error. NULL input.
#[ffi_export]
pub fn tinned_expr_is_scalar(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    match try_with_handle(h, "tinned_expr_is_scalar", "ExprHandle", |eh| {
        let expr = eh.as_ref();
        Ok(expr.is_scalar())
    }) {
        Ok(v) => v,
        Err(e) => {
            set_out_err(out_err, e);
            false
        },
    }
}

// Hash key as a string; free with `tinned_string_free`. NULL on error.
#[ffi_export]
pub fn tinned_expr_hash_key(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    match try_with_handle(h, "tinned_expr_hash_key", "ExprHandle", |eh| {
        let expr = eh.as_ref();
        Ok(tinned_string_to_cstr(expr.hash_key()))
    }) {
        Ok(s) => Some(s),
        Err(e) => {
            set_out_err(out_err, e);
            None
        },
    }
}

// Display text; free with `tinned_string_free`. NULL on error.
#[ffi_export]
pub fn tinned_expr_display(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    match try_with_handle(h, "tinned_expr_display", "ExprHandle", |eh| {
        let expr = eh.clone_arc();
        Ok(tinned_string_to_cstr(format!("{}", expr)))
    }) {
        Ok(s) => Some(s),
        Err(e) => {
            set_out_err(out_err, e);
            None
        },
    }
}

// Serialize to JSON; free with `tinned_string_free`. NULL on error.
#[ffi_export]
pub fn tinned_expr_serialize_json(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    match try_with_handle(h, "tinned_expr_serialize_json", "ExprHandle", |eh| {
        let expr = eh.as_ref();
        serde_json::to_string(expr).map(tinned_string_to_cstr).map_err(|err| {
            generic_error("Failed to serialize expression to JSON", Some(Box::new(err)))
        })
    }) {
        Ok(s) => Some(s),
        Err(e) => {
            set_out_err(out_err, e);
            None
        },
    }
}

// Deserialize from JSON. `json` is nullable `const char*`.
// Returns NULL on error. Use `tinned_expr_free` to free the result.
#[ffi_export]
pub fn tinned_expr_deserialize_json(
    json: Option<char_p::Ref<'_>>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let s = match tinned_string_from_cstr(json) {
        Some(s) => s,
        None => {
            set_out_err(
                out_err,
                generic_error("tinned_expr_deserialize_json: NULL or non-UTF-8 json", None),
            );
            return None;
        },
    };

    match serde_json::from_str::<Arc<dyn Expr>>(&s) {
        Ok(expr) => Some(ExprBox::new(ExprHandle::new(expr))),
        Err(err) => {
            set_out_err(
                out_err,
                generic_error("Failed to deserialize expression from JSON", Some(Box::new(err))),
            );
            None
        },
    }
}
