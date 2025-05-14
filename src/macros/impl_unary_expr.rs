macro_rules! impl_unary_expr_traits {
    ($type_name:ident, $type_scalar:ident, $display_fmt:expr) => {
        impl ExprInternal for $type_name {
            impl_expr_internal_methods!($type_name);

            #[inline]
            fn find_all_key(&self) -> u32 {
                self.argument.find_all_key()
            }

            #[inline]
            fn match_for_find_all(&self, other: &Arc<dyn Expr>) -> bool {
                if let Some(expr) = downcast_from_arc::<$type_name>(other) {
                    self.argument.match_for_find_all(&expr.argument)
                } else {
                    false
                }
            }

            // Except for `TwoElecOperator`, all other implemented unary `Expr`
            // types do not hold derivative. So we use equality comparison on
            // the whole expression, i.e. we do not override the method
            // `match_for_replace_all()` of the trait `ExprInternal`.
        }

        #[typetag::serde]
        impl Expr for $type_name {
            #[inline]
            fn hash_key(&self) -> String {
                format!("{}({})", stringify!($type_name), self.argument.hash_key())
            }

            impl_unary_expr_common_methods!(
                $type_name,
                argument,
                $type_scalar,
                false,
                |_this, arg| Self::new(arg)
            );

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

macro_rules! impl_unary_expr_common_methods {
    (
        $type_name:ident,
        $arg_field:ident,
        $type_scalar:ident,
        $has_derivative:tt,
        $build_expr:expr
    ) => {
        #[inline]
        fn as_any(&self) -> &dyn std::any::Any {
            self
        }

        impl_unary_expr_common_methods!(@unary_expr_is_scalar $arg_field, $type_scalar);

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
        fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
            if self.match_for_find_all(s) {
                BTreeMap::from([(self.find_all_key(), HashSet::from([self.clone_expr()]))])
            } else {
                self.$arg_field.find_all(s)
            }
        }

        #[inline]
        fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
            if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                return impl_unary_expr_common_methods!(
                    @build_zero_expr
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

        #[inline]
        fn replace(
            &self,
            map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        ) -> Result<Arc<dyn Expr>, TinnedError> {
            if let Some((_, value)) = map.iter().find(|(key, _)| self.eq_expr(key.as_ref())) {
                return Ok(value.clone());
            }

            impl_unary_expr_arg_operation!(
                self,
                $arg_field,
                |arg: &Arc<dyn Expr>| arg.replace(map),
                concat!(stringify!($type_name), "::replace() failed"),
                $build_expr
            )
        }

        #[inline]
        fn replace_all(
            &self,
            map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        ) -> Result<Arc<dyn Expr>, TinnedError> {
            if let Some((_, value)) = map.iter().find(|(key, _)| self.match_for_replace_all(key)) {
                return impl_unary_expr_common_methods!(
                    @unary_expr_replace_self
                    self,
                    value,
                    $has_derivative
                );
            }

            impl_unary_expr_arg_operation!(
                self,
                $arg_field,
                |arg: &Arc<dyn Expr>| arg.replace_all(map),
                concat!(stringify!($type_name), "::replace_all() failed"),
                $build_expr
            )
        }
    };

    (@unary_expr_is_scalar $_arg_field:ident, True) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            true
        }
    };

    (@unary_expr_is_scalar $_arg_field:ident, False) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            false
        }
    };

    (@unary_expr_is_scalar $arg_field:ident, Argument) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            self.$arg_field.is_scalar()
        }
    };

    (@unary_expr_replace_self $self:ident, $value:ident, true) => {
        if $self.derivative.is_empty() {
            Ok($value.clone())
        } else {
            differentiate_expr($value, &$self.derivative)
        }
    };

    (@unary_expr_replace_self $_self:ident, $value:ident, false) => {
        Ok($value.clone())
    };

    (@build_zero_expr $_argument:expr, True) => { impl_zero_expr!(true) };

    (@build_zero_expr $_argument:expr, False) => { impl_zero_expr!(false) };

    (@build_zero_expr $argument:expr, Argument) => {
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
