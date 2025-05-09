macro_rules! impl_unary_expr_traits {
    ($type_name:ident, $is_scalar:expr, $display_fmt:expr) => {
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

            impl_unary_expr_eq_shallow!($type_name, argument);

            #[inline]
            fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "{self}")
            }

            #[inline]
            fn clean_temporum(
                &self,
                freq_tol: Option<NumberTolerance>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                impl_unary_expr_arg_operation!(
                    self,
                    argument,
                    self.argument.clean_temporum(freq_tol),
                    concat!(stringify!($type_name), "clean_temporum() failed for argument"),
                    |arg| Self::new(arg)
                )
            }

            fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
                let diff_arg = self.argument.differentiate(s).map_err(|e| {
                    generic_expression_error(
                        concat!(stringify!($type_name), "differentiate() failed for argument"),
                        self,
                        Some(Box::new(e)),
                    )
                })?;

                Self::new(diff_arg)
            }

            #[inline]
            fn eliminate(
                &self,
                parameter: &Arc<dyn Expr>,
                perturbations: &[Arc<Perturbation>],
                min_order: u32,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                impl_unary_expr_arg_operation!(
                    self,
                    argument,
                    self.argument.eliminate(parameter, perturbations, min_order),
                    concat!(stringify!($type_name), "eliminate() failed for argument"),
                    |arg| Self::new(arg)
                )
            }

            impl_unary_expr_exist_any!(argument);
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

macro_rules! impl_unary_expr_eq_shallow {
    ($type_name:ident, $arg_field:ident) => {
        #[inline]
        fn eq_shallow(&self, other: &dyn Expr) -> bool {
            if let Some(expr) = downcast_from_ref::<$type_name>(other) {
                self.$arg_field.eq_shallow(expr.$arg_field.as_ref())
            } else {
                false
            }
        }
    };
}

macro_rules! impl_unary_expr_arg_operation {
    ($self:ident, $arg_field:ident, $arg_operation:expr, $message:expr, $build_expr:expr) => {{
        let new_arg = $arg_operation
            .map_err(|e| generic_expression_error($message, $self, Some(Box::new(e))))?;

        if &new_arg == &$self.$arg_field {
            Ok($self.clone_expr())
        } else {
            ($build_expr)(new_arg)
        }
    }};
}

macro_rules! impl_unary_expr_exist_any {
    ($arg_field:ident) => {
        #[inline]
        fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
            set.iter().any(|expr| self.eq_expr(expr.as_ref())) || self.$arg_field.exist_any(set)
        }
    };
}
