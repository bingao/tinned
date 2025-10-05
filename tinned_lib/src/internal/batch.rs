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

/// Sorts a list of expressions by grouping them by `group_key` and sorting
/// within each group by `hash_key()`.
///
/// Example:
/// let sorted_by_type = sort_expressions_grouped_by(exprs, |e| e.type_name());
/// let sorted_by_order = sort_expressions_grouped_by(exprs, |e| e.total_order());
#[inline]
pub(crate) fn sort_expressions_grouped_by<K: Ord + Copy, F: Fn(&Arc<dyn Expr>) -> K>(
    exprs: &[Arc<dyn Expr>],
    group_key: F,
) -> Vec<Arc<dyn Expr>> {
    let mut grouped: BTreeMap<K, Vec<Arc<dyn Expr>>> = BTreeMap::new();

    for expr in exprs {
        grouped.entry(group_key(expr)).or_default().push(expr.clone());
    }

    // Now flatten: within each group, sort by `hash_key()`
    let mut sorted = Vec::with_capacity(exprs.len());
    for mut group in grouped.into_values() {
        group.sort_by_key(|e| e.hash_key());
        sorted.extend(group);
    }

    sorted
}
