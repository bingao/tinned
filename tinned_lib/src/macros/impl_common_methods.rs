macro_rules! impl_expr_internal_methods {
    ($type_name:ident, $has_derivative:tt) => {
        #[inline]
        fn eq_expr(&self, other: &dyn $crate::core::Expr) -> bool {
            if let Some(expr) = $crate::public::downcast_from_ref::<$type_name>(other) {
                self == expr
            } else {
                false
            }
        }

        #[inline]
        fn fmt_expr(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
            ::std::write!(f, "{self}")
        }

        impl_expr_internal_methods!(@impl_replace_expr_self $type_name, $has_derivative);
    };

    (@impl_replace_expr_self $type_name:ident, true) => {
        #[inline]
        fn replace_expr_self(
            &self,
            expr: expr_arc_ref_ty!(),
            replacement: expr_arc_ty!(),
        ) -> expr_result_ty!() {
            if self.derivative.is_empty() {
                Ok(replacement)
            } else {
                let op = $crate::public::downcast_from_arc::<$type_name>(expr).ok_or_else(|| {
                    $crate::public::unreachable_error(
                        concat!("Expected ", stringify!($type_name)),
                        expr,
                        None,
                    )
                })?;
                $crate::public::differentiate_expr(
                    &replacement,
                    &self.derivative.complement(&op.derivative),
                )
            }
        }
    };

    (@impl_replace_expr_self $type_name:ident, false) => {};
}

macro_rules! impl_expr_common_methods {
    ($is_scalar:tt) => {
        #[inline]
        fn as_any(&self) -> &dyn ::std::any::Any {
            self
        }

        impl_is_scalar!($is_scalar);

        #[inline]
        fn clone_expr(&self) -> expr_arc_ty!() {
            ::std::sync::Arc::new(self.clone())
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
        Ok($crate::expressions::Number::zero())
    };
    (false) => {
        Ok($crate::expressions::ZeroOperator::new())
    };
    ($is_scalar:expr) => {
        if $is_scalar {
            Ok($crate::expressions::Number::zero())
        } else {
            Ok($crate::expressions::ZeroOperator::new())
        }
    };
}
