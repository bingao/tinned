use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::core::{Expr, TinnedError};
use tinned::perturbations::{PertMultichain, Perturbation};
use tinned::public::expression_error;

use crate::c_support::{tinned_string_to_cstr, try_with_handle};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::{
    PertMultichainBox, PertMultichainHandle, PerturbationBox, PerturbationHandle,
};

#[inline]
fn invalid_type_err(caller: &'static str, expr: &Arc<dyn Expr>) -> TinnedError {
    let msg: &'static str = Box::leak(
        format!("Invalid expression passed to {caller}; expected {}", expr.type_name())
            .into_boxed_str(),
    );
    expression_error(msg, expr, None)
}

// Downcast helper for value-returning closures.
// - Validates handle (NULL -> sets `out_err`, returns `None`)
// - Downcasts to `Target` (mismatch -> sets `out_err`, returns `None`)
// - Runs `f(&Target) -> R` and returns `Some(R)`
// - On any error -> sets `out_err`, returns `None`
#[inline]
pub(crate) fn with_downcast_val<Target: 'static, R: Copy>(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    caller: &'static str,
    f: impl FnOnce(&Target) -> R,
) -> Option<R> {
    try_with_handle(h, caller, "ExprHandle", |eh| {
        let expr = eh.as_ref();
        if let Some(t) = expr.as_any().downcast_ref::<Target>() {
            Ok(f(t))
        } else {
            let expr_arc = eh.clone_arc();
            Err(invalid_type_err(caller, &expr_arc))
        }
    })
    .map_or_else(
        |e| {
            tinned_error_new(out_err, e);
            None
        },
        Some,
    )
}

// Downcast helper for string-returning closures.
// - Validates handle (NULL -> sets `out_err`, returns `None`)
// - Downcasts to `Target` (mismatch -> sets `out_err`, returns `None`)
// - Runs `f(&Target) -> String` and returns `Some(char_p::Box)`
// - On any error -> sets `out_err`, returns `None`
#[inline]
pub(crate) fn with_downcast_cstr<Target: 'static>(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    caller: &'static str,
    f: impl FnOnce(&Target) -> String,
) -> Option<char_p::Box> {
    try_with_handle(h, caller, "ExprHandle", |eh| {
        let expr = eh.as_ref();
        if let Some(t) = expr.as_any().downcast_ref::<Target>() {
            Ok(tinned_string_to_cstr(f(t)))
        } else {
            let expr_arc = eh.clone_arc();
            Err(invalid_type_err(caller, &expr_arc))
        }
    })
    .map_or_else(
        |e| {
            tinned_error_new(out_err, e);
            None
        },
        Some,
    )
}

// Downcast helper for closures that produce an expression (`Arc<dyn Expr>`)
// and return a boxed handle to C.
// - Same validation/downcast as above
// - Runs `f(&Target) -> Result<Arc<dyn Expr>, TinnedError>`
// - Boxes as `ExprBox` on success
// - Sets `out_err` and returns None on failure
#[inline]
pub(crate) fn with_downcast_expr<Target: 'static>(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    caller: &'static str,
    f: impl FnOnce(&Target) -> Result<Arc<dyn Expr>, TinnedError>,
) -> Option<ExprBox> {
    try_with_handle(h, caller, "ExprHandle", |eh| {
        let expr = eh.as_ref();
        if let Some(t) = expr.as_any().downcast_ref::<Target>() {
            f(t).map(|arc| ExprBox::new(ExprHandle::new(arc)))
        } else {
            let expr_arc = eh.clone_arc();
            Err(invalid_type_err(caller, &expr_arc))
        }
    })
    .map_or_else(
        |e| {
            tinned_error_new(out_err, e);
            None
        },
        Some,
    )
}

// Downcast helper for closures that produce a perturbation (`Arc<Perturbation>`)
// and return a boxed handle to C.
// - Same validation/downcast as above
// - Runs `f(&Target) -> Result<Arc<Perturbation>, TinnedError>`
// - Boxes as `PerturbationBox` on success
// - Sets `out_err` and returns None on failure
#[inline]
pub(crate) fn with_downcast_pert<Target: 'static>(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    caller: &'static str,
    f: impl FnOnce(&Target) -> Result<Arc<Perturbation>, TinnedError>,
) -> Option<PerturbationBox> {
    try_with_handle(h, caller, "ExprHandle", |eh| {
        let expr = eh.as_ref();
        if let Some(t) = expr.as_any().downcast_ref::<Target>() {
            f(t).map(|arc| PerturbationBox::new(PerturbationHandle::new(arc)))
        } else {
            let expr_arc = eh.clone_arc();
            Err(invalid_type_err(caller, &expr_arc))
        }
    })
    .map_or_else(
        |e| {
            tinned_error_new(out_err, e);
            None
        },
        Some,
    )
}

// Downcast helper for perturbation multichain closures.
// - Validates handle (NULL -> sets `out_err`, returns `None`)
// - Downcasts to `Target` (mismatch -> sets `out_err`, returns `None`)
// - Runs `f(&Target) -> &PertMultichain` and returns `Some(PertMultichain)`
// - On any error -> sets `out_err`, returns `None`
#[inline]
pub(crate) fn with_downcast_pert_multichain<Target: 'static>(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    caller: &'static str,
    f: impl FnOnce(&Target) -> Result<Arc<PertMultichain>, TinnedError>,
) -> Option<PertMultichainBox> {
    try_with_handle(h, caller, "ExprHandle", |eh| {
        let expr = eh.as_ref();
        if let Some(t) = expr.as_any().downcast_ref::<Target>() {
            f(t).map(|arc| PertMultichainBox::new(PertMultichainHandle::new(arc)))
        } else {
            let expr_arc = eh.clone_arc();
            Err(invalid_type_err(caller, &expr_arc))
        }
    })
    .map_or_else(
        |e| {
            tinned_error_new(out_err, e);
            None
        },
        Some,
    )
}
