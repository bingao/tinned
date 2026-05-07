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

/// Free a vector of `repr_c::Vec<ExprBox>`.
/// Dropping the Vec drops each ExprBox, which decrements Arc counts.
#[ffi_export]
pub fn tinned_expr_vec_free(_v: repr_c::Vec<ExprBox>) {
    // Intentionally empty. Taking by value and returning lets _v drop here.
}

/// Borrowed slice of handles
#[repr(C)]
#[derive_ReprC]
pub struct ExprSlice {
    pub ptr: *const *const ExprHandle,
    pub len: usize,
}

/// Turn an `ExprSlice` into `Vec<Arc<dyn Expr>>`, or set `out_err` and return `None`.
#[inline]
pub fn expr_vec_from_slice(
    slice: &ExprSlice,
    caller: &'static str,
) -> Result<Vec<Arc<dyn Expr>>, TinnedError> {
    try_vec_from_slice(slice.ptr, slice.len, caller, "ExprHandle", |h: &ExprHandle| h.clone_arc())
}

// Build a HashSet<Arc<dyn Expr>> or BTreeSet<Arc<dyn Expr>> from an ExprSlice.
#[inline]
pub fn expr_set_from_slice<S>(slice: &ExprSlice, caller: &'static str) -> Result<S, TinnedError>
where
    S: Default + Extend<Arc<dyn Expr>>,
{
    try_set_from_slice(slice.ptr, slice.len, caller, "ExprHandle", |h| h.clone_arc())
}

// Build a HashMap<Arc<dyn Expr>, Arc<dyn Expr>> from parallel ExprSlices.
#[inline]
pub fn expr_map_from_slices(
    keys: &ExprSlice,
    values: &ExprSlice,
    caller: &'static str,
) -> Result<HashMap<Arc<dyn Expr>, Arc<dyn Expr>>, TinnedError> {
    try_map_from_slices(
        keys.ptr,
        keys.len,
        values.ptr,
        values.len,
        caller,
        "ExprHandle(key)",
        "ExprHandle(value)",
        |h| h.clone_arc(),
        |h| h.clone_arc(),
    )
}

// Used for `find_all` method of `Expr`
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

/// Borrowed slice of ExprSlice values.
///
/// C layout:
/// - ptr points to an array of ExprSlice
/// - len is the number of ExprSlice elements
#[repr(C)]
#[derive_ReprC]
pub struct ExprSetSlice {
    ptr: *const ExprSlice,
    len: usize,
}

#[inline]
pub fn expr_sets_from_slice(
    slice: &ExprSetSlice,
    caller: &'static str,
) -> Result<Vec<HashSet<Arc<dyn Expr>>>, TinnedError> {
    if slice.ptr.is_null() {
        if slice.len == 0 {
            return Ok(Vec::new());
        }

        return Err(generic_error(
            format!("{caller}: ExprSetSlice.ptr is NULL but len is {}", slice.len),
            None,
        ));
    }

    let slices = unsafe { std::slice::from_raw_parts(slice.ptr, slice.len) };

    slices
        .iter()
        .enumerate()
        .map(|(idx, expr_slice)| {
            expr_set_from_slice::<HashSet<Arc<dyn Expr>>>(expr_slice, caller).map_err(|e| {
                generic_error(format!("{caller}: failed to read set at index {idx}: {e}"), None)
            })
        })
        .collect()
}
