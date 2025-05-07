macro_rules! impl_unary_expr_traits {
    ($type_name:ident, $is_scalar:expr, $build_zero_expr:ident, $display_fmt:expr) => {
        #[typetag::serde]
        impl Expr for $type_name {
            #[inline]
            fn as_any(&self) -> &dyn std::any::Any {
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
            fn clone_expr(&self) -> Self {
                self.clone()
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

            #[inline]
            fn clean_temporum(
                &self,
                num_tol: Option<NumberTolerance>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                let new_arg = self.argument.clean_temporum(num_tol)?;

                if is_zero_expr(&new_arg) {
                    return Ok($build_zero_expr());
                }

                if new_arg == self.argument {
                    Ok(Arc::new(self.clone_expr()))
                } else {
                    Self::new(new_arg)
                }
            }

            fn differentiate(
                &self,
                s: &Arc<crate::perturbations::Perturbation>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                let diff_arg = self.argument.differentiate(s).map_err(|e| {
                    generic_expression_error(
                        "differentiate() on argument failed",
                        self,
                        Some(Box::new(e)),
                    )
                })?;

                Self::new(diff_arg)
            }
        }

        impl PartialEq for $type_name {
            fn eq(&self, other: &Self) -> bool {
                &self.argument == &other.argument
            }
        }

        impl Eq for $type_name {}

        impl std::fmt::Display for $type_name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, $display_fmt, arg = self.argument)
            }
        }
    };
}
