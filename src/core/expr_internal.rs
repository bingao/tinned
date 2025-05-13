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

        // Provides the key to be used for the function `find_all()`.
        #[inline]
        fn find_all_key(&self) -> u32 {
            0
        }

        // Compares equality of two expressions for the method `find_all()`.
        // This comparison will ignore derivatives of these two expression, and
        // derivatives of their fields may also be ignored. See the
        // implementation of different concrete expression types.
        #[inline]
        fn match_for_find_all(&self, other: &Arc<dyn crate::core::expr::Expr>) -> bool {
            self.eq_expr(other.as_ref())
        }

        // Compares equality of two expressions for the method `replace_all()`.
        // This comparison will ignore derivatives of these two expression, and
        // derivatives of their fields may also be ignored. See the
        // implementation of different concrete expression types.
        #[inline]
        fn match_for_replace_all(&self, other: &Arc<dyn crate::core::expr::Expr>) -> bool {
            self.eq_expr(other.as_ref())
        }
    }
}
