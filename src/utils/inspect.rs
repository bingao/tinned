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
pub fn is_expr_type<T: Any>(expr: &Arc<dyn Expr>) -> bool {
    expr.as_ref().as_any().downcast_ref::<T>().is_some()
}
