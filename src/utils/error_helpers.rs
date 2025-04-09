use std::sync::Arc;

use crate::core::{Expr, TinnedError};

#[inline]
pub fn invalid_expression_error(message: &'static str, expr: &Arc<dyn Expr>) -> TinnedError {
    TinnedError::InvalidExpression { message, expression: format!("{}", expr) }
}

#[inline]
pub fn unreachable_error(message: &'static str, expr: &Arc<dyn Expr>) -> TinnedError {
    TinnedError::Unreachable { message, expression: format!("{}", expr) }
}
