macro_rules! impl_unary_expr_traits {
    ($type_name:ident, $type_scalar:ident, $display_fmt:expr) => {
        impl $crate::core::ExprInternal for $type_name {
            impl_unary_expr_internal_methods!(
                $type_name,
                $type_scalar,
                argument,
                false,
                |_this, arg| Self::new(arg)
            );

            #[inline]
            fn hash_key(&self) -> ::std::string::String {
                ::std::format!("{}({})", stringify!($type_name), self.argument.hash_key())
            }

            #[inline]
            fn total_order(&self) -> u32 {
                self.argument.total_order()
            }

            #[inline]
            fn deep_eq_superchains(&self, other: &expr_arc_ty!()) -> bool {
                if let Some(expr) = $crate::public::downcast_from_arc::<$type_name>(other) {
                    self.argument.deep_eq_superchains(&expr.argument)
                } else {
                    false
                }
            }

            // Except for `AoTwoElecMatrix`, all other implemented unary `Expr`
            // types do not hold derivative. So we use equality comparison on
            // the whole expression, i.e. we do not override the method
            // `eq_by_superchains()` of the trait `ExprInternal`.
        }

        #[::typetag::serde]
        impl $crate::core::Expr for $type_name {
            impl_unary_expr_common_methods!($type_name, $type_scalar, argument, |_this, arg| {
                Self::new(arg)
            });

            #[inline]
            fn apply_zero_rules(
                &self,
                freq_tol: ::std::option::Option<$crate::public::NumberTolerance>,
            ) -> expr_result_ty!() {
                impl_unary_expr_arg_operation!(
                    self,
                    $type_scalar,
                    argument,
                    |arg: expr_arc_ref_ty!()| arg.apply_zero_rules(freq_tol),
                    concat!(stringify!($type_name), "::apply_zero_rules() failed"),
                    |_this, arg| Self::new(arg)
                )
            }

            #[inline]
            fn differentiate(&self, s: &pert_arc_ty!()) -> expr_result_ty!() {
                let diff_arg = self.argument.differentiate(s).map_err(|e| {
                    $crate::public::generic_expression_error(
                        concat!(stringify!($type_name), "::differentiate() failed for argument"),
                        self,
                        Some(::std::boxed::Box::new(e)),
                    )
                })?;

                Self::new(diff_arg)
            }
        }

        impl ::std::cmp::PartialEq for $type_name {
            fn eq(&self, other: &Self) -> bool {
                &self.argument == &other.argument
            }
        }

        impl ::std::cmp::Eq for $type_name {}

        impl ::std::fmt::Display for $type_name {
            fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
                ::std::write!(f, $display_fmt, arg = self.argument)
            }
        }
    };
}

macro_rules! impl_unary_expr_internal_methods {
    ($type_name:ident, $type_scalar:ident, $arg_field:ident, $has_derivative:tt, $build_expr:expr) => {
        impl_expr_internal_methods!($type_name, $has_derivative);

        #[inline]
        fn replace_expr_children(
            &self,
            map: &expr_map_ty!(),
            include_derivatives: bool,
        ) -> expr_result_ty!() {
            impl_unary_expr_arg_operation!(
                self,
                $type_scalar,
                $arg_field,
                |arg: expr_arc_ref_ty!()| arg.replace(map, include_derivatives),
                concat!(stringify!($type_name), "::replace_expr_children() failed"),
                $build_expr
            )
        }
    };
}

macro_rules! impl_unary_expr_common_methods {
    ($type_name:ident, $type_scalar:ident, $arg_field:ident, $build_expr:expr) => {
        #[inline]
        fn as_any(&self) -> &dyn ::std::any::Any {
            self
        }

        impl_unary_expr_common_methods!(@unary_is_scalar $arg_field, $type_scalar);

        #[inline]
        fn clone_expr(&self) -> expr_arc_ty!() {
            ::std::sync::Arc::new(self.clone())
        }

        #[inline]
        fn eliminate(
            &self,
            parameter: &expr_arc_ty!(),
            perturbations: &[pert_arc_ty!()],
            min_order: u32,
        ) -> expr_result_ty!() {
            impl_unary_expr_arg_operation!(
                self,
                $type_scalar,
                $arg_field,
                |arg: expr_arc_ref_ty!()| arg.eliminate(parameter, perturbations, min_order),
                concat!(stringify!($type_name), "::eliminate() failed"),
                $build_expr
            )
        }

        #[inline]
        fn exist_any(&self, set: &expr_set_ty!(), include_derivatives: bool) -> bool {
            self.match_self_any(set, include_derivatives)
                || self.$arg_field.exist_any(set, include_derivatives)
        }

        #[inline]
        fn find_superchains(
            &self,
            s: &expr_arc_ty!(),
        ) -> expr_differentiation_map_ty!() {
            if self.deep_eq_superchains(s) {
                ::std::collections::BTreeMap::from([(
                    self.total_order(),
                    ::std::collections::HashSet::from([self.clone_expr()]),
                )])
            } else {
                self.$arg_field.find_superchains(s)
            }
        }

        #[inline]
        fn remove(&self, set: &expr_set_ty!()) -> expr_result_ty!() {
            if self.match_self_any(set, false) {
                return impl_unary_zero_expr!(self.$arg_field, $type_scalar);
            }

            impl_unary_expr_arg_operation!(
                self,
                $type_scalar,
                $arg_field,
                |arg: expr_arc_ref_ty!()| arg.remove(set),
                concat!(stringify!($type_name), "::remove() failed"),
                $build_expr
            )
        }

        #[inline]
        fn retain(
            &self,
            set: &expr_set_ty!(),
            include_derivatives: bool,
        ) -> expr_result_ty!() {
            if self.match_self_any(set, include_derivatives) {
                return Ok(self.clone_expr());
            }

            impl_unary_expr_arg_operation!(
                self,
                $type_scalar,
                $arg_field,
                |arg: expr_arc_ref_ty!()| arg.retain(set, include_derivatives),
                concat!(stringify!($type_name), "::retain() failed"),
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
}

macro_rules! impl_unary_zero_expr {
    ($_argument:expr, True) => {
        impl_zero_expr!(true)
    };

    ($_argument:expr, False) => {
        impl_zero_expr!(false)
    };

    ($argument:expr, Argument) => {
        impl_zero_expr!($argument.is_scalar())
    };
}

macro_rules! impl_unary_expr_arg_operation {
    (
        $self:ident,
        $type_scalar:ident,
        $arg_field:ident,
        $arg_operation:expr,
        $message:expr,
        $build_expr:expr
    ) => {{
        let new_arg = ($arg_operation)(&$self.$arg_field).map_err(|e| {
            $crate::public::generic_expression_error(
                concat!($message, " for ", stringify!($arg_field)),
                $self,
                Some(::std::boxed::Box::new(e)),
            )
        })?;

        if $crate::public::is_zero_expr(&new_arg, None) {
            impl_unary_zero_expr!($self.$arg_field, $type_scalar)
        } else if &new_arg == &$self.$arg_field {
            Ok($self.clone_expr())
        } else {
            ($build_expr)($self, new_arg)
        }
    }};
}
