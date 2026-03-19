use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::{Debug, Display, Formatter, Result as FmtResult};
use std::hash::{Hash, Hasher};
//use log::warn;
use std::sync::Arc;

use crate::core::TinnedError;
use crate::core::expr_internal::sealed::ExprInternal;

// Base Expression trait
#[typetag::serde]
pub trait Expr: Debug + Send + Sync + ExprInternal {
    // Exposes the trait object as `dyn Any` to enable runtime downcasting of
    // trait objects to concrete types.
    fn as_any(&self) -> &dyn std::any::Any;

    // Returns the concrete type name of an expression.
    #[inline]
    fn type_name(&self) -> &'static str {
        std::any::type_name::<Self>()
    }

    // Returns hash value of the expression.
    #[inline]
    fn hash_value(&self) -> u64 {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.hash_key().hash(&mut hasher);
        hasher.finish()
    }

    // Returns if an expression is scalar.
    fn is_scalar(&self) -> bool;

    // Make a clone of an expression.
    fn clone_expr(&self) -> Arc<dyn Expr>;

    // Performs a conditional canonicalization to zero, such as setting
    // `TemporumOperator` and unperturbed `TemporumOverlap` to zero, and
    // undifferentiated perturbing operators to zero.
    #[inline]
    fn apply_zero_rules(
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

    // Checks if any expression in `set` exists in the current expression.
    #[inline]
    fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
        set.iter().any(|expr| self.eq_expr(expr.as_ref()))
    }

    // Finds a given expression `s` and all its higher-order "differentiated"
    // ones in the current expression. Here, two terms need to be clarified:
    //
    // (1) If the expression `s` has derivative, for example, `s`^{a}. Then its
    //     higher-order derivatives exclude lower-order and unrelated objects,
    //     such as `s`, `s`^{b}, or `s`^{c}, where b and c are different from
    //     a. In other words, matched expressions have derivatives that are
    //     superchains of the derivative of `s`.
    //
    // (2) The quoted term "differentiated" means derivatives may not be the
    //     mathematical derivative of `s`. For example, for `s` being the type
    //     of `TwoElecOperator`, its derivative is an `MatrixAdd` of
    //     `TwoElecOperator` objects with fields of (un)differentiated electron
    //     repulsion integrals (ERI) and one-electron spin-orbital density
    //     matrix, which is difficult to find. Instead, we return all
    //     `TwoElecOperator` objects, with (un)differentiated ERI and density
    //     matrix fields.
    //
    // To summarize, this method returns objects that match `s` according to
    // the method `deep_eq_superchains()`.
    #[inline]
    fn find_superchains(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
        if self.deep_eq_superchains(s) {
            BTreeMap::from([(self.total_order(), HashSet::from([self.clone_expr()]))])
        } else {
            BTreeMap::new()
        }
    }

    // Removes all expressions in `set` from the current expression.
    fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError>;

    // If the parameter `exact_equality` is `true`, the method replaces
    // expressions (keys of `map`) with corresponding values of `map` in the
    // current expression.
    //
    // If the parameter `exact_equality` is `false`, the method replaces
    // expressions (keys of `map`) and their higher-order "derivatives" with
    // corresponding values of `map` and their derivatives in the concrete
    // expression. Here, "higher-order" is the same as that of method
    // `find_superchains()`. The meaning of "derivatives" is taken care by
    // different concrete expression types. One requirement is that
    // `replace_superchains()` should not return same results for two different
    // `map`'s. Expressions to be replaced are determined by the method
    // `eq_by_superchains()`, which can be overriden by concrete expression
    // types.
    #[inline]
    fn replace(
        &self,
        map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        exact_equality: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let found = if exact_equality {
            map.iter()
                .find(|(key, _)| self.eq_expr(key.as_ref()))
                .map(|(_, value)| Ok(value.clone()))
        } else {
            map.iter()
                .find(|(key, _)| self.eq_by_superchains(key))
                .map(|(expr, value)| self.replace_expr_self(expr, value.clone()))
        };

        match found {
            Some(result) => result,
            None => self.replace_expr_fields(map, exact_equality),
        }
    }

    // Performs `retain_expr()` method on all expressions in `set` one by one.
    #[inline]
    fn retain(
        &self,
        set: &HashSet<Arc<dyn Expr>>,
        exact_equality: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let mut iter = set.iter();

        let first = match iter.next() {
            Some(x) => x,
            None => return Ok(self.clone_expr()),
        };

        let mut result = self.retain_expr(first, exact_equality)?;

        for expr in iter {
            if result.is_exact_zero() {
                break;
            }
            result = result.retain_expr(expr, exact_equality)?;
        }

        Ok(result)
    }

    // If the parameter `exact_equality` is `true`, the method keeps
    // sub-expressions in the current expression that contain the given `expr`,
    // while other sub-expressions are removed from the current expression.
    //
    // If the parameter `exact_equality` is `false`, sub-expressions kept
    // should contain either `expr` or its higher-order derivatives. Other
    // sub-expressions are removed even if they contain lower-order or
    // unrelated derivatives of `expr`.
    #[inline]
    fn retain_expr(
        &self,
        expr: &Arc<dyn Expr>,
        exact_equality: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let found = if exact_equality {
            self.eq_expr(expr.as_ref())
        } else {
            self.eq_by_superchains(expr)
        };

        if found {
            return Ok(self.clone_expr());
        }

        self.retain_expr_fields(expr, exact_equality)
    }
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
        // hash_key() is a string and is stable + collision-free, so we use it
        // here instead of hash_value() which is not guaranteed to be unique
        // (though hash collisions rarely happen).
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
