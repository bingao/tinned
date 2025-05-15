pub(crate) mod sealed {
    use std::fmt;
    use std::sync::Arc;

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
        // used as the key for the function `find_all()` and for sorting a list
        // of expressions. Note that "differentiation" is not mathematically
        // strict. For example, it is the differenitation only on electron
        // repulsion integrals (ERI) for `TwoElecOperator`. See implementation
        // of concrete expression types.
        #[inline]
        fn total_order(&self) -> u32 {
            0
        }

        // Compares equality of two expressions for the method `find_all()`.
        // This comparison will ignore derivatives of these two expressions, and
        // derivatives of their fields may also be ignored. See the
        // implementation of different concrete expression types.
        #[inline]
        fn match_for_find_all(&self, other: &Arc<dyn crate::core::expr::Expr>) -> bool {
            self.eq_expr(other.as_ref())
        }

        // Compares equality of two expressions for the method `replace_all()`.
        // This comparison will ignore derivatives of these two expressions.
        // For unambiguous replacement, we may require equality comparison of
        // their fields as well. See the implementation of different concrete
        // expression types.
        //
        // One requirement is that `match_for_replace_all()` should not return
        // true for two different `self`'s by only ignoring their derivatives.
        #[inline]
        fn match_for_replace_all(&self, other: &Arc<dyn crate::core::expr::Expr>) -> bool {
            self.eq_expr(other.as_ref())
        }
    }
}
