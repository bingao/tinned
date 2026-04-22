use std::any::Any;
use std::sync::Arc;

use crate::core::Expr;

// Helper function to downcast a trait object to a specific type.
// Returns `Some(&T)` if successful, otherwise `None`.
#[inline]
pub fn downcast_from_arc<T: Any>(expr: &Arc<dyn Expr>) -> Option<&T> {
    expr.as_ref().as_any().downcast_ref::<T>()
}

#[inline]
pub fn downcast_from_ref<T: Any>(expr: &dyn Expr) -> Option<&T> {
    expr.as_any().downcast_ref::<T>()
}

#[inline]
pub fn downcast_arc_expr<T>(expr: Arc<dyn Expr>) -> Option<Arc<T>>
where
    T: Expr + 'static,
{
    if !expr.as_ref().as_any().is::<T>() {
        return None;
    }

    // The runtime type check above guarantees that the allocation behind the
    // `Arc<dyn Expr>` is actually a `T`, so rebuilding the same `Arc` as
    // `Arc<T>` is valid.
    let raw = Arc::into_raw(expr) as *const T;

    unsafe { Some(Arc::from_raw(raw)) }
}

#[inline]
pub fn is_expr_type<T: Any>(expr: &Arc<dyn Expr>) -> bool {
    expr.as_ref().as_any().downcast_ref::<T>().is_some()
}
