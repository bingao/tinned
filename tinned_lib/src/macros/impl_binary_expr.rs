macro_rules! impl_binary_expr_internal_methods {
    (
        $type_name:ident,
        $scalar_rule:ident,
        $first_argument:ident,
        $second_argument:ident,
        $has_derivative:tt,
        $build_expr:expr
    ) => {
        impl_expr_internal_methods!($type_name, $has_derivative);

        #[inline]
        fn replace_one_in_children(
            &self,
            expr: expr_arc_ref_ty!(),
            replacement: expr_arc_ty!(),
            include_derivatives: bool,
        ) -> expr_result_ty!() {
            $crate::internal::transform_binary_any_zero(
                self,
                &self.$first_argument,
                &self.$second_argument,
                |arg: expr_arc_ref_ty!()| {
                    arg.replace_one(expr, replacement.clone(), include_derivatives)
                },
                concat!(stringify!($type_name), "::replace_one_in_children() failed"),
                |first_argument, second_argument| {
                    ($build_expr)(self, first_argument, second_argument)
                },
                || impl_binary_expr_zero!(self, $scalar_rule),
            )
        }

        #[inline]
        fn replace_all_in_children(
            &self,
            map: &expr_map_ty!(),
            include_derivatives: bool,
        ) -> expr_result_ty!() {
            $crate::internal::transform_binary_any_zero(
                self,
                &self.$first_argument,
                &self.$second_argument,
                |arg: expr_arc_ref_ty!()| arg.replace_all(map, include_derivatives),
                concat!(stringify!($type_name), "::replace_all_in_children() failed"),
                |first_argument, second_argument| {
                    ($build_expr)(self, first_argument, second_argument)
                },
                || impl_binary_expr_zero!(self, $scalar_rule),
            )
        }
    };
}

macro_rules! impl_binary_expr_common_methods {
    (
        $type_name:ident,
        $scalar_rule:ident,
        $first_argument:ident,
        $second_argument:ident,
        $build_expr:expr,
        $with_substitute_zero_perturbations:tt
    ) => {
        #[inline]
        fn as_any(&self) -> &dyn ::std::any::Any {
            self
        }

        #[inline]
        fn clone_expr(&self) -> expr_arc_ty!() {
            $crate::internal::intern_expr(::std::sync::Arc::new(self.clone()))
        }

        impl_binary_expr_is_scalar!($scalar_rule);

        #[inline]
        fn has_unperturbed_term(&self) -> bool {
            self.$first_argument.has_unperturbed_term()
                && self.$second_argument.has_unperturbed_term()
        }

        impl_binary_expr_common_methods!(
            @binary_expr_substitute_zero_perturbations
            $type_name,
            $scalar_rule,
            $first_argument,
            $second_argument,
            $build_expr,
            $with_substitute_zero_perturbations
        );

        #[inline]
        fn eliminate(
            &self,
            parameter: expr_arc_ty!(),
            perturbations: &[pert_arc_ty!()],
            min_order: u32,
        ) -> expr_result_ty!() {
            $crate::internal::transform_binary_any_zero(
                self,
                &self.$first_argument,
                &self.$second_argument,
                |arg: expr_arc_ref_ty!()| arg.eliminate(parameter.clone(), perturbations, min_order),
                concat!(stringify!($type_name), "::eliminate() failed"),
                |first_argument, second_argument| ($build_expr)(self, first_argument, second_argument),
                || impl_binary_expr_zero!(self, $scalar_rule),
            )
        }

        #[inline]
        fn find_all(&self, s: expr_arc_ref_ty!()) -> expr_differentiation_map_ty!() {
            if self.deep_eq_superchains(s) {
                ::std::collections::BTreeMap::from([(
                    self.total_order(),
                    ::std::collections::HashSet::from([self.clone_expr()]),
                )])
            } else {
                let mut result = self.$first_argument.find_all(s);

                for (order, subset) in self.$second_argument.find_all(s) {
                    result.entry(order).or_default().extend(subset);
                }

                result
            }
        }

        #[inline]
        fn match_one(&self, s: expr_arc_ref_ty!(), include_derivatives: bool) -> bool {
            self.match_one_self(s, include_derivatives)
                || self.$first_argument.match_one(s, include_derivatives)
                || self.$second_argument.match_one(s, include_derivatives)
        }

        #[inline]
        fn match_any(&self, set: &expr_set_ty!(), include_derivatives: bool) -> bool {
            self.match_any_self(set, include_derivatives)
                || self.$first_argument.match_any(set, include_derivatives)
                || self.$second_argument.match_any(set, include_derivatives)
        }

        #[inline]
        fn remove_one(&self, s: expr_arc_ref_ty!()) -> expr_result_ty!() {
            if self.match_one_self(s, false) {
                return impl_binary_expr_zero!(self, $scalar_rule);
            }

            $crate::internal::transform_binary_any_zero(
                self,
                &self.$first_argument,
                &self.$second_argument,
                |arg: expr_arc_ref_ty!()| arg.remove_one(s),
                concat!(stringify!($type_name), "::remove_one() failed"),
                |first_argument, second_argument| ($build_expr)(self, first_argument, second_argument),
                || impl_binary_expr_zero!(self, $scalar_rule),
            )
        }

        #[inline]
        fn remove_all(&self, set: &expr_set_ty!()) -> expr_result_ty!() {
            if self.match_any_self(set, false) {
                return impl_binary_expr_zero!(self, $scalar_rule);
            }

            $crate::internal::transform_binary_any_zero(
                self,
                &self.$first_argument,
                &self.$second_argument,
                |arg: expr_arc_ref_ty!()| arg.remove_all(set),
                concat!(stringify!($type_name), "::remove_all() failed"),
                |first_argument, second_argument| ($build_expr)(self, first_argument, second_argument),
                || impl_binary_expr_zero!(self, $scalar_rule),
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

            $crate::internal::transform_binary_all_zero(
                self,
                &self.$first_argument,
                &self.$second_argument,
                |arg: expr_arc_ref_ty!()| arg.retain_one(s, include_derivatives),
                concat!(stringify!($type_name), "::retain_one() failed"),
                |first_argument, second_argument| ($build_expr)(self, first_argument, second_argument),
                || impl_binary_expr_zero!(self, $scalar_rule),
            )
        }

        #[inline]
        fn retain_any(&self, set: &expr_set_ty!(), include_derivatives: bool) -> expr_result_ty!() {
            if self.match_any_self(set, include_derivatives) {
                return Ok(self.clone_expr());
            }

            $crate::internal::transform_binary_all_zero(
                self,
                &self.$first_argument,
                &self.$second_argument,
                |arg: expr_arc_ref_ty!()| arg.retain_any(set, include_derivatives),
                concat!(stringify!($type_name), "::retain_any() failed"),
                |first_argument, second_argument| ($build_expr)(self, first_argument, second_argument),
                || impl_binary_expr_zero!(self, $scalar_rule),
            )
        }
    };

    (@binary_expr_substitute_zero_perturbations
        $type_name:ident,
        $scalar_rule:ident,
        $first_argument:ident,
        $second_argument:ident,
        $build_expr:expr,
        true
    ) => {
        #[inline]
        fn substitute_zero_perturbations(
            &self,
            freq_tol: ::std::option::Option<$crate::public::NumberTolerance>,
        ) -> expr_result_ty!() {
            $crate::internal::transform_binary_any_zero(
                self,
                &self.$first_argument,
                &self.$second_argument,
                |arg: expr_arc_ref_ty!()| arg.substitute_zero_perturbations(freq_tol.clone()),
                concat!(stringify!($type_name), "::substitute_zero_perturbations() failed"),
                |first_argument, second_argument| ($build_expr)(self, first_argument, second_argument),
                || impl_binary_expr_zero!(self, $scalar_rule),
            )
        }
    };

    (@binary_expr_substitute_zero_perturbations
        $type_name:ident,
        $scalar_rule:ident,
        $first_argument:ident,
        $second_argument:ident,
        $build_expr:expr,
        false
    ) => {};
}

macro_rules! impl_binary_expr_zero {
    ($_self:ident, true) => {
        Ok($crate::expressions::Number::zero())
    };

    ($_self:ident, false) => {
        Ok($crate::expressions::ZeroOperator::new())
    };

    ($self:ident, FROM_SELF) => {
        if $self.is_scalar {
            Ok($crate::expressions::Number::zero())
        } else {
            Ok($crate::expressions::ZeroOperator::new())
        }
    };
}

macro_rules! impl_binary_expr_is_scalar {
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

    (FROM_SELF) => {
        #[inline]
        fn is_scalar(&self) -> bool {
            self.is_scalar
        }
    };
}
