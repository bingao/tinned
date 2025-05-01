use std::collections::BTreeMap;
use std::sync::Arc;

use crate::core::Expr;
use crate::perturbations::Perturbation;

pub(crate) fn multi_expression_format(exprs: &[Arc<dyn Expr>], delimiter: &str) -> String {
    exprs
        .iter()
        .map(|e| format!("{}", e)) // calls Display for each expression
        .collect::<Vec<_>>()
        .join(delimiter)
}

pub(crate) fn multi_expression_hash(exprs: &[Arc<dyn Expr>], delimiter: &str) -> String {
    exprs
        .iter()
        .map(|e| e.hash_key()) // calls the hash_key method for each expression
        .collect::<Vec<_>>()
        .join(delimiter)
}

pub(crate) fn multi_perturbation_format(perts: &[Arc<Perturbation>], delimiter: &str) -> String {
    perts.iter().map(|p| format!("{}", p)).collect::<Vec<_>>().join(delimiter)
}

pub(crate) fn multi_perturbation_hash(perts: &[Arc<Perturbation>], delimiter: &str) -> String {
    perts.iter().map(|p| p.hash_key()).collect::<Vec<_>>().join(delimiter)
}

/// Sorts a list of expressions by grouping them by `type_id`
/// and sorting within each group by `fast_hash()`.
#[inline]
pub(crate) fn sort_multi_expressions(exprs: &[Arc<dyn Expr>]) -> Vec<Arc<dyn Expr>> {
    // Group `exprs` by type names
    let mut grouped: BTreeMap<&'static str, Vec<Arc<dyn Expr>>> = BTreeMap::new();

    for expr in exprs {
        grouped.entry(expr.type_name()).or_default().push(expr.clone());
    }

    // Now flatten: within each group, sort by `fast_hash()`
    let mut sorted = Vec::new();
    for mut group in grouped.into_values() {
        group.sort_by_key(|e| e.fast_hash());
        sorted.extend(group);
    }

    sorted
}
