use std::any::Any;
use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter, Result as FmtResult};
use std::hash::{Hash, Hasher};
use log::warn;
use std::sync::Arc;

use crate::core::TinnedError;

// Base Expression Trait
pub trait Expr: Debug + Send + Sync {
    fn as_any(&self) -> &dyn Any;

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

    fn differentiate(
        &self,
        s: &crate::perturbations::Perturbation,
    ) -> Result<Arc<dyn Expr>, TinnedError>;
}

impl Hash for dyn Expr {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.hash_key().hash(state);
    }
}

impl PartialEq for dyn Expr {
    fn eq(&self, other: &Self) -> bool {
        // Fast path: differing hash keys => not equal
        if self.hash_key() != other.hash_key() {
            return false;
        }
        self.eq_expr(other)
    }
}

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
        Display::fmt(&*self, f)
    }
}
