macro_rules! impl_expr_internal_methods {
    ($type_name:ident) => {
        #[inline]
        fn clone_expr(&self) -> Arc<dyn Expr> {
            Arc::new(self.clone())
        }

        #[inline]
        fn eq_expr(&self, other: &dyn Expr) -> bool {
            if let Some(expr) = downcast_from_ref::<$type_name>(other) {
                self == expr
            } else {
                false
            }
        }

        #[inline]
        fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "{self}")
        }
    };
}

macro_rules! impl_expr_common_methods {
    ($is_scalar:tt) => {
        #[inline]
        fn as_any(&self) -> &dyn std::any::Any {
            self
        }

        #[inline]
        fn is_scalar(&self) -> bool {
            $is_scalar
        }
    };
}

macro_rules! impl_zero_expr {
    (true) => {
        Ok(Number::zero())
    };
    (false) => {
        Ok(ZeroOperator::new())
    };
}
