pub(crate) mod sealed {
    use std::collections::HashMap;
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
        // used as the key for the function `find_all()` and for
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
        // method `find_all()`.
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

        // Replaces the expression with `replacement`, or derivative of
        // `replacement`. `expr` is "equal to" `self` according to the method
        // `eq_by_superchains()`. So, the derivative on `replacement` can be
        // figured out by taking `complement` of the derivative on `expr`.
        #[inline]
        fn apply_replacement(
            &self,
            _expr: &Arc<dyn crate::core::expr::Expr>,
            replacement: Arc<dyn crate::core::expr::Expr>,
        ) -> Result<Arc<dyn crate::core::expr::Expr>, TinnedError> {
            Ok(replacement)
        }

        // Performs `replace_one()` method on child subexpression(s) if there
        // exists By default, there is no child subexpression and we simply
        // return the clone of the expression.
        fn replace_one_in_children(
            &self,
            _expr: &Arc<dyn crate::core::expr::Expr>,
            _replacement: Arc<dyn crate::core::expr::Expr>,
            _include_derivatives: bool,
        ) -> Result<Arc<dyn crate::core::expr::Expr>, TinnedError>;

        // Performs `replace_all()` method on child subexpression(s) if there
        // exists By default, there is no child subexpression and we simply
        // return the clone of the expression.
        fn replace_all_in_children(
            &self,
            _map: &HashMap<Arc<dyn crate::core::expr::Expr>, Arc<dyn crate::core::expr::Expr>>,
            _include_derivatives: bool,
        ) -> Result<Arc<dyn crate::core::expr::Expr>, TinnedError>;

        // Returns if the expression is exactly zero
        #[inline]
        fn is_exact_zero(&self) -> bool {
            false
        }
    }
}
