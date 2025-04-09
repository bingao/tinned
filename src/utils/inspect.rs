use std::sync::Arc;

use crate::core::Expr;
use crate::expressions::{Number, ZeroOperator};

// Helper function to downcast a trait object to a specific type.
// Returns `Some(&T)` if successful, otherwise `None`.
#[inline]
pub fn downcast_from_arc<T: std::any::Any>(expr: &Arc<dyn Expr>) -> Option<&T> {
    expr.as_ref().as_any().downcast_ref::<T>()
}

#[inline]
pub fn downcast_from_ref<T: std::any::Any>(expr: &dyn Expr) -> Option<&T> {
    expr.as_any().downcast_ref::<T>()
}

#[inline]
pub fn is_zero_expr(expr: &Arc<dyn Expr>) -> bool {
    downcast_from_arc::<Number>(expr).map_or(false, |n| n.is_zero()) || expr.is::<ZeroOperator>()
}

#[inline]
pub fn is_one_expr(expr: &Arc<dyn Expr>) -> bool {
    downcast_from_arc::<Number>(expr).map_or(false, |n| n.is_one())
}
