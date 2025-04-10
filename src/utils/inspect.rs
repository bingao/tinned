use std::any::Any;
use std::sync::Arc;

use crate::core::Expr;
use crate::expressions::{Number, ZeroOperator};

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

#[inline]
pub fn is_zero_expr(expr: &Arc<dyn Expr>) -> bool {
    expr.as_any().downcast_ref::<ZeroOperator>().is_some()
        || downcast_from_arc::<Number>(expr).map_or(false, |n| n.is_zero())
}

#[inline]
pub fn is_one_expr(expr: &Arc<dyn Expr>) -> bool {
    downcast_from_arc::<Number>(expr).map_or(false, |n| n.is_one())
}
