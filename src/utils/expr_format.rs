use std::sync::Arc;

use crate::core::Expr;

pub fn join_exprs_for_display(exprs: &[Arc<dyn Expr>], delimiter: &str) -> String {
    exprs
        .iter()
        .map(|e| format!("{}", e)) // calls Display for each expression
        .collect::<Vec<_>>()
        .join(delimiter)
}

pub fn join_exprs_for_hash(exprs: &[Arc<dyn Expr>], delimiter: &str) -> String {
    exprs
        .iter()
        .map(|e| e.hash_key()) // calls the hash_key method for each expression
        .collect::<Vec<_>>()
        .join(delimiter)
}
