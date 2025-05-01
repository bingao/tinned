use std::sync::Arc;

use crate::core::Expr;
use crate::perturbations::Perturbation;

// An interner is a global cache that:
// * Stores previously created expressions to reuse them instead of re-allocating.
// * Reduces memory usage by avoiding duplicate storage of the same expression.
macro_rules! define_interner {
    ($fn_name:ident, $ty:ty, $map_name:ident) => {
        lazy_static::lazy_static! {
            static ref $map_name: Arc<dashmap::DashMap<String, Arc<$ty>>> =
                Arc::new(dashmap::DashMap::new());
        }

        #[inline]
        pub(crate) fn $fn_name(obj: Arc<$ty>) -> Arc<$ty> {
            let key = obj.hash_key();

            match $map_name.entry(key) {
                dashmap::mapref::entry::Entry::Occupied(e) => Arc::clone(e.get()),
                dashmap::mapref::entry::Entry::Vacant(e) => {
                    e.insert(Arc::clone(&obj));
                    obj
                },
            }
        }
    };
}

define_interner!(intern_expr, dyn Expr, EXPR_INTERNER);
define_interner!(intern_pert, Perturbation, PERT_INTERNER);
