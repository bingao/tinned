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

    // Returns name of a concrete expression type.
    #[inline]
    fn type_name(&self) -> &'static str {
        std::any::type_name::<Self>()
    }

    // Returns hash key of a concrete expression type.
    fn hash_key(&self) -> String;

    // Precomputes hashes for faster sorting.
    #[inline]
    fn fast_hash(&self) -> u64 {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.hash_key().hash(&mut hasher);
        hasher.finish()
    }

    // Returns if a concrete expression type is scalar.
    fn is_scalar(&self) -> bool;

    // Compares equality for concrete expression types.
    fn eq_expr(&self, other: &dyn Expr) -> bool;

    // Compares equality for concrete expression types but ignores derivatives.
    fn eq_shallow(&self, other: &Arc<dyn Expr>) -> bool {
        self.eq_expr(other.as_ref())
    }

    // Make a clone of a concrete expression type.
    fn clone_expr(&self) -> Arc<dyn Expr>;

    // Formats a concrete expression type.
    fn fmt_expr(&self, f: &mut Formatter) -> FmtResult;

    // Cleans `TemporumOperator` and unperturbed `TemporumOverlap` objects.
    #[inline]
    fn clean_temporum(
        &self,
        _freq_tol: Option<crate::public::NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(self.clone_expr())
    }

    // Differentiates with respect to a `Perturbation`.
    fn differentiate(
        &self,
        s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError>;

    // Eliminates a given response `parameter`'s derivatives from the
    // expression. Maximum order of derivatives to be eliminated is the length
    // of `perturbations`, and minimum order is specified by `min_order`. For
    // wave function parameters, it should be greater than the floor function
    // of the half length of `perturbations`, and for multipliers, it should be
    // greater than or equal to the ceiling function of the half length of
    // `perturbations` according to J. Chem. Phys. 129, 214103 (2008).
    //
    // Note that we expect that `parameter` is either `LagMultiplier` or
    // `WfnParameter`. Error or incorrect result will return if users provide
    // invalid types of `parameter`.
    #[inline]
    fn eliminate(
        &self,
        _parameter: &Arc<dyn Expr>,
        _perturbations: &[Arc<crate::perturbations::Perturbation>],
        _min_order: u32,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(self.clone_expr())
    }

    // Checks if any expression in `set` exists in the concrete expression.
    #[inline]
    fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
        set.iter().any(|expr| self.eq_expr(expr.as_ref()))
    }

    //    // Finds a given expression `s` and all its differentiated ones in the
    //    // concrete expression.
    //    fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>>;

    //    // Removes given expressions in `set` from the concrete expression.
    //    fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError>;

    //    // Replaces given expressions (keys of `map`) and their derivatives with
    //    // corresponding values of `map` and their derivatives in the concrete
    //    // expression.
    //    fn replace(&self, map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>) -> Arc<dyn Expr>;
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
