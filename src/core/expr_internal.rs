pub(crate) mod sealed {
    use std::collections::{HashMap, HashSet};
    use std::fmt;
    use std::sync::Arc;

    use crate::core::TinnedError;

    pub trait ExprInternal {
        // Make a clone of an expression.
        fn clone_expr(&self) -> Arc<dyn crate::core::expr::Expr>;

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
        // on electron repulsion integrals (ERI) for `TwoElecOperator`. See
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

        // Performs `replace()` method on the expression's field(s) if it has
        // any. By default, there is no fields in the expression and we simply
        // return its clone.
        fn replace_expr_fields(
            &self,
            _map: &HashMap<Arc<dyn crate::core::expr::Expr>, Arc<dyn crate::core::expr::Expr>>,
            _exact_equality: bool,
        ) -> Result<Arc<dyn crate::core::expr::Expr>, TinnedError> {
            Ok(self.clone_expr())
        }

        // Performs `retain()` method on the expression's field(s) if it has any.
        fn retain_expr_fields(
            &self,
            set: &HashSet<Arc<dyn crate::core::expr::Expr>>,
            exact_equality: bool,
        ) -> Result<Arc<dyn crate::core::expr::Expr>, TinnedError>;
    }
}
