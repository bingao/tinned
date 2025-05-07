use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::{Debug, Display, Formatter, Result as FmtResult};
use std::hash::{Hash, Hasher};
//use log::warn;
use std::sync::Arc;

use typetag;

use crate::core::TinnedError;

// Base Expression trait
#[typetag::serde]
pub trait Expr: Debug + Send + Sync {
    //
    fn as_any(&self) -> &dyn std::any::Any;

    // Returns name of a concrete expression type
    #[inline]
    fn type_name(&self) -> &'static str {
        std::any::type_name::<Self>()
    }

    // Returns hash key of a concrete expression type
    fn hash_key(&self) -> String;

    // Precomputes hashes for faster sorting
    #[inline]
    fn fast_hash(&self) -> u64 {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.hash_key().hash(&mut hasher);
        hasher.finish()
    }

    // Returns if a concrete expression type is scalar
    fn is_scalar(&self) -> bool;

    // Compares equality for concrete expression types
    fn eq_expr(&self, other: &dyn Expr) -> bool;

    // Make a clone of a concrete expression type
    fn clone_expr(&self) -> Self
    where
        Self: Sized + Clone;

    // Formats a concrete expression type
    fn fmt_expr(&self, f: &mut Formatter) -> FmtResult;

    // Cleans `TemporumOperator` and unperturbed `TemporumOverlap` objects
    #[inline]
    fn clean_temporum(
        &self,
        _num_tol: Option<crate::public::NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(Arc::new(self.clone_expr()))
    }

    // Differentiate with respect to a `Perturbation`
    fn differentiate(
        &self,
        s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError>;

    // Helper function to eliminate a given response `parameter`'s derivatives
    // from `x`. Maximum order of derivatives to be eliminated is the length of
    // `perturbations`, and minimum order is specified by `min_order`. For wave
    // function parameters, it should be greater than the floor function of the
    // half length of `perturbations`, and for multipliers, it should be greater
    // than or equal to the ceiling function of the half length of
    // `perturbations` according to J. Chem. Phys. 129, 214103 (2008).
    fn eliminate(
        &self,
        parameter: &Arc<dyn Expr>,
        perturbations: &[Arc<Perturbation>],
        min_order: u32,
    ) -> Result<Arc<dyn Expr>, TinnedError>;

    //
    fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool;

    // Find a given `s` and all its differentiated ones
    fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>>;

    // Helper function to remove given `symbols` from `x`
    fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError>;

    // Helper function to replace Tinned objects and their derivatives with
    // SymEngine `Basic` symbols and corresponding derivatives
    fn replace(&self, map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>) -> Arc<dyn Expr>;
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
