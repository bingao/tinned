macro_rules! impl_unary_expr_traits {
    ($type_name:ident, $type_scalar:ident, $display_fmt:expr) => {
        impl ExprInternal for $type_name {
            impl_unary_expr_internal_methods!($type_name, argument, false, |_this, arg| Self::new(
                arg
            ));

            #[inline]
            fn hash_key(&self) -> String {
                format!("{}({})", stringify!($type_name), self.argument.hash_key())
            }

            #[inline]
            fn total_order(&self) -> u32 {
                self.argument.total_order()
            }

            #[inline]
            fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
                if let Some(expr) = downcast_from_arc::<$type_name>(other) {
                    self.argument.deep_eq_superchains(&expr.argument)
                } else {
                    false
                }
            }

            // Except for `TwoElecOperator`, all other implemented unary `Expr`
            // types do not hold derivative. So we use equality comparison on
            // the whole expression, i.e. we do not override the method
            // `eq_by_superchains()` of the trait `ExprInternal`.
        }

        #[typetag::serde]
        impl Expr for $type_name {
            impl_unary_expr_common_methods!($type_name, argument, $type_scalar, |_this, arg| {
                Self::new(arg)
            });

            #[inline]
            fn clean_temporum(
                &self,
                freq_tol: Option<NumberTolerance>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                impl_unary_expr_arg_operation!(
                    self,
                    argument,
                    |arg: &Arc<dyn Expr>| arg.clean_temporum(freq_tol),
                    concat!(stringify!($type_name), "::clean_temporum() failed"),
                    |_this, arg| Self::new(arg)
                )
            }

            #[inline]
            fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
                let diff_arg = self.argument.differentiate(s).map_err(|e| {
                    generic_expression_error(
                        concat!(stringify!($type_name), "::differentiate() failed for argument"),
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

macro_rules! impl_unary_expr_internal_methods {
    ($type_name:ident, $arg_field:ident, $has_derivative:tt, $build_expr:expr) => {
        impl_expr_internal_methods!($type_name, $has_derivative);

        #[inline]
        fn replace_expr_fields(
            &self,
            map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
            exact_equality: bool,
        ) -> Result<Arc<dyn Expr>, TinnedError> {
            impl_unary_expr_arg_operation!(
                self,
                $arg_field,
                |arg: &Arc<dyn Expr>| arg.replace(map, exact_equality),
                concat!(stringify!($type_name), "::replace_expr_fields() failed"),
                $build_expr
            )
        }

        #[inline]
        fn retain_expr_fields(
            &self,
            set: &HashSet<Arc<dyn Expr>>,
            exact_equality: bool,
        ) -> Result<Arc<dyn Expr>, TinnedError> {
            impl_unary_expr_arg_operation!(
                self,
                $arg_field,
                |arg: &Arc<dyn Expr>| arg.retain(set, exact_equality),
                concat!(stringify!($type_name), "::retain_expr_fields() failed"),
                $build_expr
            )
        }
    };
}

macro_rules! impl_unary_expr_common_methods {
    ($type_name:ident, $arg_field:ident, $type_scalar:ident, $build_expr:expr) => {
        #[inline]
        fn as_any(&self) -> &dyn std::any::Any {
            self
        }

        impl_unary_expr_common_methods!(@unary_is_scalar $arg_field, $type_scalar);

        #[inline]
        fn clone_expr(&self) -> Arc<dyn Expr> {
            Arc::new(self.clone())
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
                $arg_field,
                |arg: &Arc<dyn Expr>| arg.eliminate(parameter, perturbations, min_order),
                concat!(stringify!($type_name), "::eliminate() failed"),
                $build_expr
            )
        }

        #[inline]
        fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
            set.iter().any(|expr| self.eq_expr(expr.as_ref())) || self.$arg_field.exist_any(set)
        }

        #[inline]
        fn find_superchains(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
            if self.deep_eq_superchains(s) {
                BTreeMap::from([(self.total_order(), HashSet::from([self.clone_expr()]))])
            } else {
                self.$arg_field.find_superchains(s)
            }
        }

        #[inline]
        fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
            if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                return impl_unary_expr_common_methods!(
                    @unary_build_zero_expr
                    self.$arg_field,
                    $type_scalar
                );
            }

            impl_unary_expr_arg_operation!(
                self,
                $arg_field,
                |arg: &Arc<dyn Expr>| arg.remove(set),
                concat!(stringify!($type_name), "::remove() failed"),
                $build_expr
            )
        }
    };

    (@unary_is_scalar $_arg_field:ident, True) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            true
        }
    };

    (@unary_is_scalar $_arg_field:ident, False) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            false
        }
    };

    (@unary_is_scalar $arg_field:ident, Argument) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            self.$arg_field.is_scalar()
        }
    };

    (@unary_build_zero_expr $_argument:expr, True) => { impl_zero_expr!(true) };

    (@unary_build_zero_expr $_argument:expr, False) => { impl_zero_expr!(false) };

    (@unary_build_zero_expr $argument:expr, Argument) => {
        if $argument.is_scalar() {
            impl_zero_expr!(true)
        } else {
            impl_zero_expr!(false)
        }
    };
}

macro_rules! impl_unary_expr_arg_operation {
    (
        $self:ident,
        $arg_field:ident,
        $arg_operation:expr,
        $message:expr,
        $build_expr:expr
    ) => {{
        let new_arg = ($arg_operation)(&$self.$arg_field).map_err(|e| {
            generic_expression_error(
                concat!($message, " for ", stringify!($arg_field)),
                $self,
                Some(Box::new(e)),
            )
        })?;

        if &new_arg == &$self.$arg_field {
            Ok($self.clone_expr())
        } else {
            ($build_expr)($self, new_arg)
        }
    }};
}
