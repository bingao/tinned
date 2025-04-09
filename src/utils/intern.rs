use dashmap::DashMap;
use std::sync::Arc;

use crate::core::Expr;

// An interner is a global cache that:
// * Stores previously created expressions to reuse them instead of re-allocating.
// * Reduces memory usage by avoiding duplicate storage of the same expression.
lazy_static::lazy_static! {
    static ref INTERNER: Arc<DashMap<String, Arc<dyn Expr>>> = Arc::new(DashMap::new());
}

#[inline]
pub fn intern(expr: Arc<dyn Expr>) -> Arc<dyn Expr> {
    let key = expr.hash_key();

    if let Some(cached) = INTERNER.get(&key) {
        return Arc::clone(&cached);
    }

    INTERNER.insert(key.clone(), Arc::clone(&expr));
    expr
}
