use std::sync::Arc;

use crate::core::Expr;
use crate::internal::intern_expr;

pub fn expr_to_json(expr: &Arc<dyn Expr>) -> Result<String, serde_json::Error> {
    serde_json::to_string(expr)
}

pub fn expr_from_json(s: &str) -> Result<Arc<dyn Expr>, serde_json::Error> {
    let expr: Arc<dyn Expr> = serde_json::from_str(s)?;
    Ok(intern_expr(expr))
}
