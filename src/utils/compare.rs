use std::sync::Arc;

use crate::core::Expr;
use crate::expressions::{Number, ZeroOperator};
use crate::utils::{NumberTolerance, downcast_from_arc};

#[inline]
pub fn is_zero_expr(expr: &Arc<dyn Expr>, num_tol: Option<NumberTolerance>) -> bool {
    expr.as_any().downcast_ref::<ZeroOperator>().is_some()
        || downcast_from_arc::<Number>(expr).map_or(false, |n| n.is_zero(num_tol))
}

#[inline]
pub fn is_one_expr(expr: &Arc<dyn Expr>, num_tol: Option<NumberTolerance>) -> bool {
    downcast_from_arc::<Number>(expr).map_or(false, |n| n.is_one(num_tol))
}
