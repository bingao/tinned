use std::sync::Arc;

use tinned::core::Expr;

use crate::c_support::{vec_arc_from_ptrs, with_box_or_err};
use crate::core::TinnedErrorBox;

// Opaque wrappers that C holds as pointers
pub struct ExprBox {
    inner: Arc<dyn Expr>,
}

impl ExprBox {
    #[inline]
    pub(crate) fn new(expr: Arc<dyn Expr>) -> Self {
        Self {
            inner: expr,
        }
    }

    // Clones the inner `Arc` (refcount +1) and return it.
    #[inline]
    pub(crate) fn arc_clone(&self) -> Arc<dyn Expr> {
        Arc::clone(&self.inner)
    }

    // Borrows the inner trait object (Rust-only; never expose to C).
    // #[inline]
    // pub(crate) fn as_ref(&self) -> &(dyn Expr + 'static) {
    //     self.inner.as_ref()
    // }

    // Borrow the owning `Arc` (needed when downstream APIs require `&Arc<dyn Expr>`).
    #[inline]
    pub(crate) fn as_arc(&self) -> &Arc<dyn Expr> {
        &self.inner
    }
}

// Allocates from an existing `Arc<dyn Expr>`
#[inline]
pub(crate) fn expr_box_from(expr: Arc<dyn Expr>) -> *mut ExprBox {
    Box::into_raw(Box::new(ExprBox::new(expr)))
}

// Borrows `&Arc<dyn Expr>` safely from a raw handle. Returns `None` if `h` is `NULL`.
#[inline]
pub(crate) fn with_expr_or_err<R>(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
    caller: &'static str,
    f: impl FnOnce(&Arc<dyn Expr>) -> R,
) -> Option<R> {
    with_box_or_err::<ExprBox, Arc<dyn Expr>, R>(h, out_err, caller, "ExprBox", ExprBox::as_arc, f)
}

// Converts an array of `ExprBox` to `Vec<Arc<dyn Expr>>`.
// - `ptrs` must point to an array of `count` elements of type `*const ExprBox`,
//   properly aligned and alive for the duration of this call.
// - Each element may be null (we check and error out), otherwise must point to a valid `ExprBox`.
#[inline]
pub(crate) unsafe fn vec_expr_from_ptrs(
    ptrs: *const *const ExprBox,
    count: usize,
    caller: &'static str,
    out_err: *mut *mut TinnedErrorBox,
) -> Option<Vec<Arc<dyn Expr>>> {
    unsafe {
        vec_arc_from_ptrs::<ExprBox, dyn Expr>(
            ptrs,
            count,
            caller,
            "Expr",
            out_err,
            ExprBox::arc_clone,
        )
    }
}
