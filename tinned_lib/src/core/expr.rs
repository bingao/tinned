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

    // Makes a clone of an expression.
    fn clone_expr(&self) -> Arc<dyn Expr>;

    // Returns if an expression is scalar.
    fn is_scalar(&self) -> bool;

    // Returns if an expression has unperturbed term, or zeroth-order term.
    #[inline]
    fn has_unperturbed_term(&self) -> bool {
        true
    }

    // Sets all perturbations to zero, which includes, such as setting
    // `TimeEvolution` and unperturbed `BasisTimeEvolution` to zero, and
    // undifferentiated perturbing operators to zero.
    //
    // The frequency tolerance `freq_tol` is used to determine whether a
    // numerical frequency can be treated zero or not.
    #[inline]
    fn substitute_zero_perturbations(
        &self,
        _freq_tol: Option<crate::public::NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(self.clone_expr())
    }

    // Differentiates with respect to a `Perturbation`.
    fn differentiate(
        &self,
        s: Arc<crate::perturbations::Perturbation>,
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
        _parameter: Arc<dyn Expr>,
        _perturbations: &[Arc<crate::perturbations::Perturbation>],
        _min_order: u32,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(self.clone_expr())
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
    //     of `AoTwoElecMatrix`, its derivative is an `MatrixAdd` of
    //     `AoTwoElecMatrix` objects with fields of (un)differentiated electron
    //     repulsion integrals (ERI) and one-electron spin-orbital density
    //     matrix, which is difficult to find. Instead, we return all
    //     `AoTwoElecMatrix` objects, with (un)differentiated ERI and density
    //     matrix fields.
    //
    // To summarize, this method returns objects that match `s` according to
    // the method `deep_eq_superchains()`.
    #[inline]
    fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
        if self.deep_eq_superchains(s) {
            BTreeMap::from([(self.expr_order(), HashSet::from([self.clone_expr()]))])
        } else {
            BTreeMap::new()
        }
    }

    // For a given expression `s`, checks if the current expression is equal to
    // `s`, or is a derivative (including order 0) of `s` when
    // `include_derivatives` is `true`.
    //
    //FIXME: implement its FFI
    #[inline]
    fn match_one_self(&self, s: &Arc<dyn Expr>, include_derivatives: bool) -> bool {
        if include_derivatives {
            self.eq_by_superchains(s)
        } else {
            self.eq_expr(s.as_ref())
        }
    }

    // This method is similar to `match_one_self()`, but when the current
    // expression does not match the expression `s`, this method will
    // call `match_one()` on all its child expressions if exist. This
    // method can be viewed as a recursive version of matching, and should be
    // overridden by concrete expression struct's when they have child
    // expression(s).
    //
    //FIXME: implement its FFI
    #[inline]
    fn match_one(&self, s: &Arc<dyn Expr>, include_derivatives: bool) -> bool {
        self.match_one_self(s, include_derivatives)
    }

    // For a given `set` of expressions, checks if the current expression is
    // equal to, or a derivative (including order 0) of any expression in the
    // `set`. Matching derivatives is enabled when `include_derivatives` is
    // `true`.
    //
    //FIXME: implement its FFI
    #[inline]
    fn match_any_self(&self, set: &HashSet<Arc<dyn Expr>>, include_derivatives: bool) -> bool {
        if include_derivatives {
            set.iter().any(|expr| self.eq_by_superchains(expr))
        } else {
            set.iter().any(|expr| self.eq_expr(expr.as_ref()))
        }
    }

    // This method is similar to `match_any_self()`, but when the current
    // expression does not match any expression in the `set`, this method will
    // call `match_any()` on all its child expressions if exist. This
    // method can be viewed as a recursive version of matching, and should be
    // overridden by concrete expression struct's when they have child
    // expression(s).
    #[inline]
    fn match_any(&self, set: &HashSet<Arc<dyn Expr>>, include_derivatives: bool) -> bool {
        self.match_any_self(set, include_derivatives)
    }

    // Removes a given expressions `s` from the current expression.
    fn remove_one(&self, s: &Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>;

    // Removes all expressions in `set` from the current expression.
    fn remove_all(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError>;

    // If the parameter `include_derivatives` is `false`, the method replaces
    // `expr` with `replacement` in the current expression.
    //
    // If the parameter `include_derivatives` is `true`, the method replaces
    // `expr` and its higher-order "derivatives" with corresponding values of
    // `replacement` and its derivatives in the current expression. Here,
    // "higher-order" is the same as that of method `find_all()`. The meaning
    // of "derivatives" is taken care by different concrete expression types.
    // One requirement is that `replace_one()` should not return same results
    // for two different paris of `expr` and `replacement`. Expressions to be
    // replaced are determined by the method `eq_by_superchains()`, which can
    // be overriden by concrete expression types.
    //
    //FIXME: add its FFI
    #[inline]
    fn replace_one(
        &self,
        expr: &Arc<dyn Expr>,
        replacement: Arc<dyn Expr>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if include_derivatives {
            if self.eq_by_superchains(expr) {
                return self.apply_replacement(expr, replacement);
            }
        } else if self.eq_expr(expr.as_ref()) {
            return Ok(replacement);
        }

        self.replace_one_in_children(expr, replacement, include_derivatives)
    }

    // If the parameter `include_derivatives` is `false`, the method replaces
    // expressions (keys of `map`) with corresponding values of `map` in the
    // current expression.
    //
    // If the parameter `include_derivatives` is `true`, the method replaces
    // expressions (keys of `map`) and their higher-order "derivatives" with
    // corresponding values of `map` and their derivatives in the current
    // expression. Here, "higher-order" is the same as that of method
    // `find_all()`. The meaning of "derivatives" is taken care by
    // different concrete expression types. One requirement is that
    // `replace_all()` should not return same results for two different
    // `map`'s. Expressions to be replaced are determined by the method
    // `eq_by_superchains()`, which can be overriden by concrete expression
    // types.
    #[inline]
    fn replace_all(
        &self,
        map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        // Note that we cannot iterate `map` and use the method `replace_one()`
        // because we will first check if the current expression matches any
        // key of `map` and make the replacement if matching. Unless no
        // matching at all, we will work on its child expression(s).
        let found = if include_derivatives {
            map.iter()
                .find(|(key, _)| self.eq_by_superchains(key))
                .map(|(expr, value)| self.apply_replacement(expr, value.clone()))
        } else {
            map.iter()
                .find(|(key, _)| self.eq_expr(key.as_ref()))
                .map(|(_, value)| Ok(value.clone()))
        };

        match found {
            Some(result) => result,
            None => self.replace_all_in_children(map, include_derivatives),
        }
    }

    /// Retains parts of the expression that match the given expression `s`.
    ///
    /// If the current expression matches `s`, or (when
    /// `include_derivatives` is `true`) / corresponds to a higher-order
    /// derivative of `s`, it is kept unchanged.
    ///
    /// Otherwise, if the current expression has no child expressions, zero
    /// is returned. If it has child expressions, the same procedure is
    /// applied recursively to each child. Based on the results, the
    /// function may return zero, the original expression, or a modified
    /// expression, depending on how the concrete expression type combines
    /// its children.
    ///
    /// This function is composable and may be applied repeatedly, for
    /// example in higher-order residue computations.
    ///
    /// FIXME: add its FFI
    fn retain_one(
        &self,
        s: &Arc<dyn crate::core::expr::Expr>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn crate::core::expr::Expr>, TinnedError>;

    /// Retains parts of the expression that match any expression in the `set`.
    ///
    /// If the current expression matches at least one expression in the `set`,
    /// or (when `include_derivatives` is `true`) corresponds to a
    /// higher-order derivative of the expression, it is kept unchanged.
    ///
    /// Otherwise, if the current expression has no child expressions, zero
    /// is returned. If it has child expressions, the same procedure is
    /// applied recursively to each child. Based on the results, the
    /// function may return zero, the original expression, or a modified
    /// expression, depending on how the concrete expression type combines
    /// its children.
    //
    //FIXME: add its FFI
    fn retain_any(
        &self,
        set: &HashSet<Arc<dyn Expr>>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError>;

    /// Applies successive retention operations with respect to each set of
    /// expressions in `sets`.
    ///
    /// Starting from the current expression, this function repeatedly applies
    /// [`retain_any`] for each set of expressions, updating the expression
    /// at each step. Conceptually, this corresponds to extracting the
    /// component of the expression that is consistent with all retention
    /// conditions induced by elements of `sets`.
    ///
    /// When `include_derivatives` is `true`, higher-order derivatives are also
    /// considered in each retention step.
    ///
    /// If at any stage the expression becomes exactly zero, the procedure
    /// terminates early and zero is returned.
    ///
    /// This function is primarily intended for (higher-order) residue
    /// computations.
    fn retain_all(
        &self,
        sets: &[HashSet<Arc<dyn Expr>>],
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let mut result = self.clone_expr();

        for set in sets {
            result = result.retain_any(set, include_derivatives)?;

            if result.is_exact_zero() {
                return Ok(result);
            }
        }

        Ok(result)
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
