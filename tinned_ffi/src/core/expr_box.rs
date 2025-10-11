use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::core::{Expr, TinnedError};

use crate::c_support::try_from_slice;

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
    pub(crate) fn new(expr: Arc<dyn Expr>) -> Self {
        Self {
            inner: expr,
        }
    }

    #[inline]
    pub(crate) fn as_ref(&self) -> &dyn Expr {
        &*self.inner
    }

    #[inline]
    pub(crate) fn clone_arc(&self) -> Arc<dyn Expr> {
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
    try_from_slice(slice, caller, "ExprHandle", |h: &ExprHandle| h.clone_arc())
}
