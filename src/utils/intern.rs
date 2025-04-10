use crate::core::Expr;
use crate::perturbations::Perturbation;

// An interner is a global cache that:
// * Stores previously created expressions to reuse them instead of re-allocating.
// * Reduces memory usage by avoiding duplicate storage of the same expression.
macro_rules! define_interner {
    ($fn_name:ident, $ty:ty, $map_name:ident) => {
        lazy_static::lazy_static! {
            static ref $map_name: std::sync::Arc<dashmap::DashMap<String, std::sync::Arc<$ty>>> =
                std::sync::Arc::new(dashmap::DashMap::new());
        }

        #[inline]
        pub fn $fn_name(obj: std::sync::Arc<$ty>) -> std::sync::Arc<$ty> {
            use std::sync::Arc;

            let key = obj.hash_key();

            if let Some(cached) = $map_name.get(&key) {
                return Arc::clone(&cached);
            }

            $map_name.insert(key.clone(), Arc::clone(&obj));
            obj
        }
    };
}

define_interner!(intern_expr, dyn Expr, EXPR_INTERNER);
define_interner!(intern_pert, Perturbation, PERT_INTERNER);
