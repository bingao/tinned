use std::collections::BTreeMap;
use std::sync::Arc;

use crate::core::Expr;

// Join mapped strings
pub(crate) fn join_mapped<I, T, F>(items: I, delimiter: &str, map_fn: F) -> String
where
    I: IntoIterator<Item = T>,
    F: Fn(T) -> String,
{
    let mut out = String::new();

    for (i, item) in items.into_iter().enumerate() {
        if i > 0 {
            out.push_str(delimiter);
        }
        out.push_str(&map_fn(item));
    }

    out
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
