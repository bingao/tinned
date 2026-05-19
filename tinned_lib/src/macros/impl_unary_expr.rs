macro_rules! impl_unary_expr_traits {
    ($type_name:ident, $scalar_rule:ident, $display_fmt:expr) => {
        impl $crate::core::ExprInternal for $type_name {
            impl_unary_expr_internal_methods!(
                $type_name,
                $scalar_rule,
                argument,
                false,
                |_this, arg| Self::new(arg)
            );

            #[inline]
            fn hash_key(&self) -> ::std::string::String {
                ::std::format!("{}({})", stringify!($type_name), self.argument.hash_key())
            }

            #[inline]
            fn expr_order(&self) -> u32 {
                self.argument.expr_order()
            }

            #[inline]
            fn deep_eq_superchains(&self, other: expr_arc_ref_ty!()) -> bool {
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
            impl_unary_expr_common_methods!($type_name, $scalar_rule, argument, |_this, arg| {
                Self::new(arg)
            });

            #[inline]
            fn has_unperturbed_term(&self) -> bool {
                self.argument.has_unperturbed_term()
            }

            #[inline]
            fn substitute_zero_perturbations(
                &self,
                freq_tol: ::std::option::Option<$crate::public::NumberTolerance>,
            ) -> expr_result_ty!() {
                $crate::internal::transform_unary_any_zero(
                    self,
                    &self.argument,
                    |arg: expr_arc_ref_ty!()| arg.substitute_zero_perturbations(freq_tol),
                    concat!(
                        stringify!($type_name),
                        "::substitute_zero_perturbations() failed for argument"
                    ),
                    |arg| Self::new(arg),
                    || impl_unary_expr_zero!(self.argument, $scalar_rule),
                )
            }

            #[inline]
            fn differentiate(&self, s: pert_arc_ty!()) -> expr_result_ty!() {
                $crate::internal::transform_unary_any_zero(
                    self,
                    &self.argument,
                    |arg: expr_arc_ref_ty!()| arg.differentiate(s),
                    concat!(stringify!($type_name), "::differentiate() failed for argument"),
                    |arg| Self::new(arg),
                    || impl_unary_expr_zero!(self.argument, $scalar_rule),
                )
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
    ($type_name:ident, $scalar_rule:ident, $arg_field:ident, $has_derivative:tt, $build_expr:expr) => {
        impl_expr_internal_methods!($type_name, $has_derivative);

        #[inline]
        fn replace_one_in_children(
            &self,
            expr: expr_arc_ref_ty!(),
            replacement: expr_arc_ty!(),
            include_derivatives: bool,
        ) -> expr_result_ty!() {
            $crate::internal::transform_unary_any_zero(
                self,
                &self.$arg_field,
                |arg: expr_arc_ref_ty!()| arg.replace_one(expr, replacement, include_derivatives),
                concat!(
                    stringify!($type_name),
                    "::replace_one_in_children() failed for ",
                    stringify!($arg_field)
                ),
                |arg| ($build_expr)(self, arg),
                || impl_unary_expr_zero!(self.$arg_field, $scalar_rule),
            )
        }

        #[inline]
        fn replace_all_in_children(
            &self,
            map: &expr_map_ty!(),
            include_derivatives: bool,
        ) -> expr_result_ty!() {
            $crate::internal::transform_unary_any_zero(
                self,
                &self.$arg_field,
                |arg: expr_arc_ref_ty!()| arg.replace_all(map, include_derivatives),
                concat!(
                    stringify!($type_name),
                    "::replace_all_in_children() failed for ",
                    stringify!($arg_field)
                ),
                |arg| ($build_expr)(self, arg),
                || impl_unary_expr_zero!(self.$arg_field, $scalar_rule),
            )
        }
    };
}

macro_rules! impl_unary_expr_common_methods {
    ($type_name:ident, $scalar_rule:ident, $arg_field:ident, $build_expr:expr) => {
        #[inline]
        fn as_any(&self) -> &dyn ::std::any::Any {
            self
        }

        #[inline]
        fn clone_expr(&self) -> expr_arc_ty!() {
            $crate::internal::intern_expr(::std::sync::Arc::new(self.clone()))
        }

        impl_unary_expr_is_scalar!($arg_field, $scalar_rule);

        #[inline]
        fn eliminate(
            &self,
            parameter: expr_arc_ty!(),
            perturbations: &[pert_arc_ty!()],
            min_order: u32,
        ) -> expr_result_ty!() {
            $crate::internal::transform_unary_any_zero(
                self,
                &self.$arg_field,
                |arg: expr_arc_ref_ty!()| arg.eliminate(parameter, perturbations, min_order),
                concat!(
                    stringify!($type_name),
                    "::eliminate() failed for ",
                    stringify!($arg_field)
                ),
                |arg| ($build_expr)(self, arg),
                || impl_unary_expr_zero!(self.$arg_field, $scalar_rule),
            )
        }

        #[inline]
        fn find_all(&self, s: expr_arc_ref_ty!()) -> expr_differentiation_map_ty!() {
            if self.deep_eq_superchains(s) {
                ::std::collections::BTreeMap::from([(
                    self.expr_order(),
                    ::std::collections::HashSet::from([self.clone_expr()]),
                )])
            } else {
                self.$arg_field.find_all(s)
            }
        }

        #[inline]
        fn match_one(&self, s: expr_arc_ref_ty!(), include_derivatives: bool) -> bool {
            self.match_one_self(s, include_derivatives)
                || self.$arg_field.match_one(s, include_derivatives)
        }

        #[inline]
        fn match_any(&self, set: &expr_set_ty!(), include_derivatives: bool) -> bool {
            self.match_any_self(set, include_derivatives)
                || self.$arg_field.match_any(set, include_derivatives)
        }

        #[inline]
        fn remove_one(&self, s: expr_arc_ref_ty!()) -> expr_result_ty!() {
            if self.match_one_self(s, false) {
                return impl_unary_expr_zero!(self.$arg_field, $scalar_rule);
            }

            $crate::internal::transform_unary_any_zero(
                self,
                &self.$arg_field,
                |arg: expr_arc_ref_ty!()| arg.remove_one(s),
                concat!(
                    stringify!($type_name),
                    "::remove_one() failed for ",
                    stringify!($arg_field)
                ),
                |arg| ($build_expr)(self, arg),
                || impl_unary_expr_zero!(self.$arg_field, $scalar_rule),
            )
        }

        #[inline]
        fn remove_all(&self, set: &expr_set_ty!()) -> expr_result_ty!() {
            if self.match_any_self(set, false) {
                return impl_unary_expr_zero!(self.$arg_field, $scalar_rule);
            }

            $crate::internal::transform_unary_any_zero(
                self,
                &self.$arg_field,
                |arg: expr_arc_ref_ty!()| arg.remove_all(set),
                concat!(
                    stringify!($type_name),
                    "::remove_all() failed for ",
                    stringify!($arg_field)
                ),
                |arg| ($build_expr)(self, arg),
                || impl_unary_expr_zero!(self.$arg_field, $scalar_rule),
            )
        }

        #[inline]
        fn retain_one(
            &self,
            s: expr_arc_ref_ty!(),
            include_derivatives: bool,
        ) -> expr_result_ty!() {
            if self.match_one_self(s, include_derivatives) {
                return Ok(self.clone_expr());
            }

            $crate::internal::transform_unary_any_zero(
                self,
                &self.$arg_field,
                |arg: expr_arc_ref_ty!()| arg.retain_one(s, include_derivatives),
                concat!(
                    stringify!($type_name),
                    "::retain_one() failed for ",
                    stringify!($arg_field)
                ),
                |arg| ($build_expr)(self, arg),
                || impl_unary_expr_zero!(self.$arg_field, $scalar_rule),
            )
        }

        #[inline]
        fn retain_any(&self, set: &expr_set_ty!(), include_derivatives: bool) -> expr_result_ty!() {
            if self.match_any_self(set, include_derivatives) {
                return Ok(self.clone_expr());
            }

            $crate::internal::transform_unary_any_zero(
                self,
                &self.$arg_field,
                |arg: expr_arc_ref_ty!()| arg.retain_any(set, include_derivatives),
                concat!(
                    stringify!($type_name),
                    "::retain_any() failed for ",
                    stringify!($arg_field)
                ),
                |arg| ($build_expr)(self, arg),
                || impl_unary_expr_zero!(self.$arg_field, $scalar_rule),
            )
        }
    };
}

macro_rules! impl_unary_expr_zero {
    ($_arg_field:expr, true) => {
        Ok($crate::expressions::Number::zero())
    };

    ($_arg_field:expr, false) => {
        Ok($crate::expressions::ZeroOperator::new())
    };

    ($arg_field:expr, FROM_ARG) => {
        if $arg_field.is_scalar() {
            Ok($crate::expressions::Number::zero())
        } else {
            Ok($crate::expressions::ZeroOperator::new())
        }
    };
}

macro_rules! impl_unary_expr_is_scalar {
    ($_arg_field:ident, true) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            true
        }
    };

    ($_arg_field:ident, false) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            false
        }
    };

    ($arg_field:ident, FROM_ARG) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            self.$arg_field.is_scalar()
        }
    };
}
