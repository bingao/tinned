use std::sync::Arc;

use crate::core::Expr;
use crate::perturbations::Perturbation;

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

pub fn join_perts_for_display(perts: &[Arc<Perturbation>], delimiter: &str) -> String {
    perts.iter().map(|p| format!("{}", p)).collect::<Vec<_>>().join(delimiter)
}

pub fn join_perts_for_hash(perts: &[Arc<Perturbation>], delimiter: &str) -> String {
    perts.iter().map(|p| p.hash_key()).collect::<Vec<_>>().join(delimiter)
}
