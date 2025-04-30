use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter, Result as FmtResult};
use std::hash::{Hash, Hasher};
//use log::warn;
use std::sync::Arc;

use typetag;

use crate::core::TinnedError;

// Base Expression Trait
#[typetag::serde]
pub trait Expr: Debug + Send + Sync {
    fn as_any(&self) -> &dyn std::any::Any;

    #[inline]
    fn type_name(&self) -> &'static str {
        std::any::type_name::<Self>()
    }

    fn hash_key(&self) -> String;

    // Precompute Hashes for Faster Sorting
    #[inline]
    fn fast_hash(&self) -> u64 {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.hash_key().hash(&mut hasher);
        hasher.finish()
    }

    fn is_scalar(&self) -> bool;

    // Compare equality for concrete expression types
    fn eq_expr(&self, other: &dyn Expr) -> bool;

    fn fmt_expr(&self, f: &mut Formatter) -> FmtResult;

    fn differentiate(
        &self,
        s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError>;
}

impl Hash for dyn Expr {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.hash_key().hash(state);
    }
}

impl PartialEq for dyn Expr {
    fn eq(&self, other: &Self) -> bool {
        // Shortcut for exact matches
        if self.hash_key() == other.hash_key() {
            return true;
        }
        // Otherwise fall back to semantic comparison
        self.eq_expr(other)
    }
}

impl Eq for dyn Expr {}

impl PartialOrd for dyn Expr {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.hash_key().cmp(&other.hash_key()))
    }
}

impl Ord for dyn Expr {
    fn cmp(&self, other: &Self) -> Ordering {
        self.hash_key().cmp(&other.hash_key())
    }
}

impl Display for dyn Expr {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        self.fmt_expr(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    test_struct_safety!(Arc<dyn Expr>);
}
