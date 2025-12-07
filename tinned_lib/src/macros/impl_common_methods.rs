macro_rules! impl_expr_internal_methods {
    ($type_name:ident, $has_derivative:tt) => {
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

        impl_expr_internal_methods!(@impl_replace_expr_self $type_name, $has_derivative);
    };

    (@impl_replace_expr_self $type_name:ident, true) => {
        #[inline]
        fn replace_expr_self(
            &self,
            expr: &Arc<dyn Expr>,
            replacement: Arc<dyn Expr>,
        ) -> Result<Arc<dyn Expr>, TinnedError> {
            if self.derivative.is_empty() {
                Ok(replacement)
            } else {
                let op = downcast_from_arc::<$type_name>(expr).ok_or_else(|| {
                    unreachable_error(concat!("Expected ", stringify!($type_name)), expr, None)
                })?;
                differentiate_expr(&replacement, &self.derivative.complement(&op.derivative))
            }
        }
    };

    (@impl_replace_expr_self $type_name:ident, false) => { }
}

macro_rules! impl_expr_common_methods {
    ($is_scalar:tt) => {
        #[inline]
        fn as_any(&self) -> &dyn std::any::Any {
            self
        }

        impl_is_scalar!($is_scalar);

        #[inline]
        fn clone_expr(&self) -> Arc<dyn Expr> {
            Arc::new(self.clone())
        }
    };
}

macro_rules! impl_is_scalar {
    (true) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            true
        }
    };
    (false) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            false
        }
    };
    ($is_scalar:tt) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            self.$is_scalar
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
    ($is_scalar:expr) => {
        if $is_scalar {
            Ok(Number::zero())
        } else {
            Ok(ZeroOperator::new())
        }
    };
}
