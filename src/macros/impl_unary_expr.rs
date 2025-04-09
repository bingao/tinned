macro_rules! impl_unary_expr_traits {
    ($type_name:ident, $is_scalar:expr, $display_fmt:expr) => {
        impl Expr for $type_name {
            #[inline]
            fn as_any(&self) -> &dyn Any {
                self
            }

            #[inline]
            fn hash_key(&self) -> String {
                format!("{}({})", stringify!($type_name), self.argument.hash_key())
            }

            #[inline]
            fn is_scalar(&self) -> bool {
                $is_scalar
            }

            #[inline]
            fn eq_expr(&self, other: &dyn Expr) -> bool {
                if let Some(expr) = downcast_expr::<$type_name>(other) {
                    self.argument == expr.argument
                } else {
                    false
                }
            }

            fn differentiate(&self, s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError> {
                let diff_arg = self.argument.differentiate(s)?;
                Self::new(diff_arg)
            }
        }

        impl Display for $type_name {
            fn fmt(&self, f: &mut Formatter) -> FmtResult {
                write!(f, $display_fmt, arg = self.argument)
            }
        }
    };
}
