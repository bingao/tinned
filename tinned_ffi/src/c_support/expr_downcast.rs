use safer_ffi::prelude::*;
use std::collections::BTreeSet;
use std::sync::Arc;

use tinned::core::{Expr, TinnedError};
use tinned::inspect::downcast_from_ref;
use tinned::perturbations::Perturbation;
use tinned::public::expression_error;

use crate::c_support::try_with_handle;
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::{PerturbationBox, PerturbationHandle};

#[inline]
fn invalid_expr_type(caller: &'static str, expr: &Arc<dyn Expr>) -> TinnedError {
    let msg: &'static str = Box::leak(
        format!("Invalid expression passed to {caller}; expected {}", expr.type_name())
            .into_boxed_str(),
    );
    expression_error(msg, expr, None)
}

// Expr downcast helper
// - Validates handle (NULL -> sets `out_err`, returns `None`)
// - Downcasts to `Target` (mismatch -> sets `out_err`, returns `None`)
// - Runs `f(&Target) -> Result<R, TinnedError>` and returns `Some(R)`
// - On any error -> sets `out_err`, returns `None`
#[inline]
pub(crate) fn ffi_map_expr_as<Target: 'static, R>(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    caller: &'static str,
    f: impl FnOnce(&Target) -> Result<R, TinnedError>,
) -> Option<R> {
    let res = try_with_handle(h, caller, "ExprHandle", |eh| {
        let expr = eh.as_ref();
        if let Some(t) = downcast_from_ref::<Target>(expr) {
            f(t)
        } else {
            let expr_arc = eh.clone_arc();
            Err(invalid_expr_type(caller, &expr_arc))
        }
    });

    match res {
        Ok(v) => Some(v),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Expr downcast helper for for small-copy returns.
#[inline]
pub(crate) fn ffi_map_expr_as_copy<Target: 'static, R: Copy>(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    caller: &'static str,
    f: impl FnOnce(&Target) -> R,
) -> Option<R> {
    ffi_map_expr_as::<Target, R>(h, out_err, caller, |t| Ok(f(t)))
}

// Expr downcast helper for for &[Arc<dyn Expr>] returns.
#[inline]
pub(crate) fn ffi_map_expr_as_exprvec<Target: 'static>(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    caller: &'static str,
    to_slice: impl FnOnce(&Target) -> &[Arc<dyn Expr>],
) -> repr_c::Vec<ExprBox> {
    match ffi_map_expr_as::<Target, _>(h, out_err, caller, |t| {
        let slice = to_slice(t);

        // Build a standard Vec<ExprBox> first
        let mut v: Vec<ExprBox> = Vec::with_capacity(slice.len());
        for expr_arc in slice {
            v.push(ExprBox::new(ExprHandle::new(Arc::clone(expr_arc))));
        }

        // Convert to repr_c::Vec<ExprBox>
        Ok(v.into())
    }) {
        Some(v) => v,
        None => Vec::<ExprBox>::new().into(),
    }
}

// Expr downcast helper for for &BTreeSet<Arc<Perturbation>> returns.
#[inline]
pub(crate) fn ffi_map_expr_as_pertvec<Target: 'static>(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    caller: &'static str,
    to_set: impl FnOnce(&Target) -> &BTreeSet<Arc<Perturbation>>,
) -> repr_c::Vec<PerturbationBox> {
    match ffi_map_expr_as::<Target, _>(h, out_err, caller, |t| {
        let set = to_set(t);

        let mut v: Vec<PerturbationBox> = Vec::with_capacity(set.len());
        for pert_arc in set {
            v.push(PerturbationBox::new(PerturbationHandle::new(Arc::clone(pert_arc))));
        }

        Ok(v.into())
    }) {
        Some(v) => v,
        None => Vec::<PerturbationBox>::new().into(),
    }
}
