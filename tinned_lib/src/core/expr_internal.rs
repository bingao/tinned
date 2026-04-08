pub(crate) mod sealed {
    use std::collections::{HashMap, HashSet};
    use std::fmt;
    use std::sync::Arc;

    use crate::core::TinnedError;

    pub trait ExprInternal {
        // Compares equality for two expressions.
        fn eq_expr(&self, other: &dyn crate::core::expr::Expr) -> bool;

        // Formats an expression.
        fn fmt_expr(&self, f: &mut fmt::Formatter) -> fmt::Result;

        // Returns hash key of an expression.
        fn hash_key(&self) -> String;

        // Total order of "differentiation" on the expression, which can be
        // used as the key for the function `find_superchains()` and for
        // sorting a list of expressions. Note that "differentiation" is not
        // mathematically strict. For example, it is the differenitation only
        // on electron repulsion integrals (ERI) for `AoTwoElecMatrix`. See
        // implementation of concrete expression types.
        #[inline]
        fn total_order(&self) -> u32 {
            0
        }

        // This equality comparison usually requires `self`'s derivative is the
        // superchain of that of `other`. Expression fields are usually
        // compared by using the same rule, while non-expression fields usually
        // requires exact equality. See the implementation of different
        // concrete expression types. This comparison is mostly used by the
        // method `find_superchains()`.
        #[inline]
        fn deep_eq_superchains(&self, other: &Arc<dyn crate::core::expr::Expr>) -> bool {
            self.eq_expr(other.as_ref())
        }

        // This equality comparison usually requires `self`'s derivative is the
        // superchain of that of `other`, while exact equality of other fields.
        // See the implementation of different concrete expression types.
        #[inline]
        fn eq_by_superchains(&self, other: &Arc<dyn crate::core::expr::Expr>) -> bool {
            self.eq_expr(other.as_ref())
        }

        //
        #[inline]
        fn match_self_single(
            &self,
            s: &Arc<dyn crate::core::expr::Expr>,
            include_derivatives: bool,
        ) -> bool {
            if include_derivatives {
                self.eq_by_superchains(s)
            } else {
                self.eq_expr(s.as_ref())
            }
        }

        // Checks if any expression in `set` exists in the current expression
        // `self`. If the parameter `include_derivatives` is `true`, we also
        // consider derivatives (including order 0) of expressions in `set`
        // when checking existence. Different from `exist_any()`, this function
        // will not check existence for the children of `self`.
        //
        // This function will be used by `retain()` as well.
        #[inline]
        fn match_self_any(
            &self,
            set: &HashSet<Arc<dyn crate::core::expr::Expr>>,
            include_derivatives: bool,
        ) -> bool {
            if include_derivatives {
                set.iter().any(|expr| self.eq_by_superchains(expr))
            } else {
                set.iter().any(|expr| self.eq_expr(expr.as_ref()))
            }
        }

        // Replaces the expression with `replacement`, or derivative of
        // `replacement`. `expr` is "equal to" `self` according to the method
        // `eq_by_superchains()`. So, the derivative on `replacement` can be
        // figured out by taking `complement` of the derivative on `expr`.
        #[inline]
        fn replace_expr_self(
            &self,
            _expr: &Arc<dyn crate::core::expr::Expr>,
            replacement: Arc<dyn crate::core::expr::Expr>,
        ) -> Result<Arc<dyn crate::core::expr::Expr>, TinnedError> {
            Ok(replacement)
        }

        // Performs `replace()` method on child subexpression(s) if there
        // exists By default, there is no child subexpression and we simply
        // return the clone of the expression.
        fn replace_expr_children(
            &self,
            _map: &HashMap<Arc<dyn crate::core::expr::Expr>, Arc<dyn crate::core::expr::Expr>>,
            _include_derivatives: bool,
        ) -> Result<Arc<dyn crate::core::expr::Expr>, TinnedError>;

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
        fn retain_single(
            &self,
            s: &Arc<dyn crate::core::expr::Expr>,
            include_derivatives: bool,
        ) -> Result<Arc<dyn crate::core::expr::Expr>, TinnedError>;

        // Returns if the expression is exactly zero
        #[inline]
        fn is_exact_zero(&self) -> bool {
            false
        }
    }
}
