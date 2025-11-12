use safer_ffi::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use tinned::core::{Expr, TinnedError};
use tinned::public::generic_error;

use crate::c_support::{try_map_from_slices, try_set_from_slice, try_vec_from_slice};
use crate::core::{TinnedErrorBox, tinned_error_new};

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

// Free an expression (NULL-safe).
#[ffi_export]
pub fn tinned_expr_free(expr: Option<ExprBox>) {
    drop(expr);
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

// Used for `find_superchains` method of `Expr`
#[derive_ReprC]
#[repr(opaque)]
pub struct ExprSuperchainHandle {
    // Sorted derivative orders
    orders: Vec<u32>,
    // For each `orders[i]`, the expressions at that order
    order_exprs: Vec<Vec<Arc<dyn Expr>>>,
}

impl ExprSuperchainHandle {
    #[inline]
    pub fn new(orders: Vec<u32>, order_exprs: Vec<Vec<Arc<dyn Expr>>>) -> Self {
        Self {
            orders,
            order_exprs,
        }
    }
}

// C-owns this box once it's returned.
pub type ExprSuperchainBox = repr_c::Box<ExprSuperchainHandle>;

#[ffi_export]
pub fn tinned_expr_superchains_len(h: Option<&ExprSuperchainHandle>) -> usize {
    h.map(|h| h.orders.len()).unwrap_or(0)
}

#[ffi_export]
pub fn tinned_expr_superchains_order_at(
    h: Option<&ExprSuperchainHandle>,
    order_idx: usize,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> u32 {
    match h.and_then(|h| h.orders.get(order_idx)).copied() {
        Some(order) => order,
        None => {
            tinned_error_new(
                out_err,
                generic_error(format!("Index {order_idx} out of range"), None),
            );
            0
        },
    }
}

#[ffi_export]
pub fn tinned_expr_superchains_order_len(
    h: Option<&ExprSuperchainHandle>,
    order_idx: usize,
) -> usize {
    h.and_then(|h| h.order_exprs.get(order_idx)).map(|v| v.len()).unwrap_or(0)
}

#[ffi_export]
pub fn tinned_expr_superchains_order_expr_at(
    h: Option<&ExprSuperchainHandle>,
    order_idx: usize,
    expr_idx: usize,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    match h.and_then(|h| h.order_exprs.get(order_idx)).and_then(|v| v.get(expr_idx)) {
        Some(expr) => Some(ExprBox::new(ExprHandle::new(expr.clone()))),
        None => {
            tinned_error_new(
                out_err,
                generic_error(
                    format!("Index out of range (order_idx={}, expr_idx={})", order_idx, expr_idx),
                    None,
                ),
            );
            None
        },
    }
}

/// Free a superchains object (NULL-safe).
#[ffi_export]
pub fn tinned_expr_superchains_free(h: Option<ExprSuperchainBox>) {
    // Taking by value transfers ownership back to Rust; drop runs on return.
    drop(h);
}
