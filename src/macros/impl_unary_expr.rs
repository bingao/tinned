macro_rules! impl_unary_expr_traits {
    ($type_name:ident, $type_scalar:ident, $display_fmt:expr) => {
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
                |_this, arg| Self::new(arg),
                true,
                true,
            );

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
        $build_expr:expr,
        $with_eq_shallow:tt,
        $with_clean_temporum:tt,
    ) => {
        #[inline]
        fn as_any(&self) -> &dyn std::any::Any {
            self
        }

        impl_unary_expr_common_methods!(@unary_expr_is_scalar $type_scalar);

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

        impl_unary_expr_common_methods!(
            @unary_expr_eq_shallow
            $type_name,
            $arg_field,
            $with_eq_shallow
        );

        #[inline]
        fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "{self}")
        }

        #[inline]
        fn total_order(&self) -> u32 {
            self.$arg_field.total_order()
        }

        impl_unary_expr_common_methods!(
            @unary_expr_clean_temporum
            $type_name,
            $arg_field,
            $build_expr,
            $with_clean_temporum
        );

        #[inline]
        fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
            set.iter().any(|expr| self.eq_expr(expr.as_ref())) || self.$arg_field.exist_any(set)
        }

        #[inline]
        fn eliminate(
            &self,
            parameter: &Arc<dyn Expr>,
            perturbations: &[Arc<Perturbation>],
            min_order: u32,
        ) -> Result<Arc<dyn Expr>, TinnedError> {
            impl_unary_expr_common_methods!(
                @unary_expr_arg_operation
                self,
                $arg_field,
                self.$arg_field.eliminate(parameter, perturbations, min_order),
                concat!(stringify!($type_name), "eliminate() failed for argument"),
                $build_expr,
            )
        }

        #[inline]
        fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
            if self.eq_shallow(s) {
                BTreeMap::from([(self.total_order(), HashSet::from([self.clone_expr()]))])
            } else {
                self.$arg_field.find_all(s)
            }
        }
    };

    (@unary_expr_is_scalar True) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            true
        }
    };

    (@unary_expr_is_scalar False) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            false
        }
    };

    (@unary_expr_is_scalar Argument) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            self.argument.is_scalar()
        }
    };

    (@unary_expr_eq_shallow $type_name:ident, $arg_field:ident, true) => {
        #[inline]
        fn eq_shallow(&self, other: &Arc<dyn Expr>) -> bool {
            if let Some(expr) = downcast_from_arc::<$type_name>(other) {
                self.$arg_field.eq_shallow(&expr.$arg_field)
            } else {
                false
            }
        }

    };

    (@unary_expr_eq_shallow $type_name:ident, $arg_field:ident, false) => { };

    (@unary_expr_clean_temporum $type_name:ident, $arg_field:ident, $build_expr:expr, true) => {
        #[inline]
        fn clean_temporum(
            &self,
            freq_tol: Option<NumberTolerance>,
        ) -> Result<Arc<dyn Expr>, TinnedError> {
            impl_unary_expr_common_methods!(
                @unary_expr_arg_operation
                self,
                $arg_field,
                self.$arg_field.clean_temporum(freq_tol),
                concat!(stringify!($type_name), "clean_temporum() failed for argument"),
                $build_expr,
            )
        }
    };

    (@unary_expr_clean_temporum $type_name:ident, $arg_field:ident, $build_expr:expr, false) => { };

    (@unary_expr_arg_operation
        $self:ident,
        $arg_field:ident,
        $arg_operation:expr,
        $message:expr,
        $build_expr:expr,
    ) => {{
        let new_arg = $arg_operation
            .map_err(|e| generic_expression_error($message, $self, Some(Box::new(e))))?;

        if &new_arg == &$self.$arg_field {
            Ok($self.clone_expr())
        } else {
            ($build_expr)($self, new_arg)
        }
    }};
}
