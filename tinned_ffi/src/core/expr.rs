use safer_ffi::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use tinned::core::{Expr, TinnedError};

use crate::c_support::{try_map_from_slices, try_set_from_slice, try_vec_from_slice};

/// An *opaque* handle that C can only pass around
#[derive_ReprC]
#[repr(opaque)]
pub struct ExprHandle {
    inner: Arc<dyn Expr>,
}

/// Owned box by C after return
pub type ExprBox = repr_c::Box<ExprHandle>;

impl ExprHandle {
    #[inline]
    pub fn new(expr: Arc<dyn Expr>) -> Self {
        Self {
            inner: expr,
        }
    }

    #[inline]
    pub fn as_ref(&self) -> &dyn Expr {
        &*self.inner
    }

    #[inline]
    pub fn clone_arc(&self) -> Arc<dyn Expr> {
        Arc::clone(&self.inner)
    }
}

/// Borrowed slice of handles
pub type ExprSlice<'a> = c_slice::Ref<'a, *const ExprHandle>;

/// Turn an `ExprSlice` into `Vec<Arc<dyn Expr>>`, or set `out_err` and return `None`.
#[inline]
pub fn expr_vec_from_slice(
    slice: ExprSlice<'_>,
    caller: &'static str,
) -> Result<Vec<Arc<dyn Expr>>, TinnedError> {
    try_vec_from_slice(slice, caller, "ExprHandle", |h: &ExprHandle| h.clone_arc())
}

// Build a HashSet<Arc<dyn Expr>> from an ExprSlice.
#[inline]
pub fn expr_set_from_slice(
    slice: ExprSlice<'_>,
    caller: &'static str,
) -> Result<HashSet<Arc<dyn Expr>>, TinnedError> {
    try_set_from_slice::<ExprHandle, dyn Expr>(slice, caller, "ExprHandle", |h| h.clone_arc())
}

// Build a HashMap<Arc<dyn Expr>, Arc<dyn Expr>> from parallel ExprSlices.
#[inline]
pub fn expr_map_from_slices(
    keys: ExprSlice<'_>,
    values: ExprSlice<'_>,
    caller: &'static str,
) -> Result<HashMap<Arc<dyn Expr>, Arc<dyn Expr>>, TinnedError> {
    try_map_from_slices::<ExprHandle, ExprHandle, dyn Expr, dyn Expr>(
        keys,
        values,
        caller,
        "ExprHandle(key)",
        "ExprHandle(value)",
        |h| h.clone_arc(),
        |h| h.clone_arc(),
    )
}

// Free an expression (NULL-safe).
#[ffi_export]
pub fn tinned_expr_free(expr: Option<ExprBox>) {
    drop(expr);
}
