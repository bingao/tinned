use std::collections::BTreeMap;
use std::sync::Arc;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Add, MatrixAdd, MatrixMul, Mul, Number};
use crate::utils::{downcast_from_arc, invalid_expression_error};

#[inline]
pub fn negate_expr(expr: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
    if let Some(num) = downcast_from_arc::<Number>(&expr) {
        return Ok(num.negate().into());
    }

    let minus_one = Number::minus_one();

    if expr.is_scalar() {
        Mul::new(vec![minus_one, expr])
    } else {
        MatrixMul::new(vec![minus_one, expr])
    }
}

#[inline]
pub fn subtract_exprs(
    lhs: Arc<dyn Expr>,
    rhs: Arc<dyn Expr>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    match (lhs.is_scalar(), rhs.is_scalar()) {
        (true, true) => Add::new(vec![lhs, negate_expr(rhs)?]),
        (false, false) => MatrixAdd::new(vec![lhs, negate_expr(rhs)?]),
        _ => Err(invalid_expression_error(
            "subtract_exprs: lhs and rhs must both be scalar or both be non-scalar",
            &lhs,
        )),
    }
}

/// Sorts a list of expressions by grouping them by `type_id`
/// and sorting within each group by `fast_hash()`.
#[inline]
pub(crate) fn group_and_sort_terms(terms: Vec<Arc<dyn Expr>>) -> Vec<Arc<dyn Expr>> {
    // Group terms by type names
    let mut grouped: BTreeMap<&'static str, Vec<Arc<dyn Expr>>> = BTreeMap::new();

    for term in terms {
        grouped.entry(term.type_name()).or_default().push(term);
    }

    // Now flatten: within each group, sort by `fast_hash()`
    let mut sorted = Vec::new();
    for mut group in grouped.into_values() {
        group.sort_by_key(|term| term.fast_hash());
        sorted.extend(group);
    }

    sorted
}
