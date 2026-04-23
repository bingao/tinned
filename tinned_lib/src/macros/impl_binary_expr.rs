macro_rules! impl_binary_expr_internal_methods {
    (
        $type_name:ident,
        $is_scalar:tt,
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
            impl_binary_expr_arg_operation!(
                self,
                $is_scalar,
                $first_argument,
                $second_argument,
                |first_arg: expr_arc_ref_ty!()| first_arg.replace_one(
                    expr,
                    replacement.clone(),
                    include_derivatives
                ),
                |second_arg: expr_arc_ref_ty!()| second_arg.replace_one(
                    expr,
                    replacement,
                    include_derivatives
                ),
                concat!(stringify!($type_name), "::replace_one_in_children() failed"),
                $build_expr
            )
        }

        #[inline]
        fn replace_all_in_children(
            &self,
            map: &expr_map_ty!(),
            include_derivatives: bool,
        ) -> expr_result_ty!() {
            impl_binary_expr_arg_operation!(
                self,
                $is_scalar,
                $first_argument,
                $second_argument,
                |arg: expr_arc_ref_ty!()| arg.replace_all(map, include_derivatives),
                concat!(stringify!($type_name), "::replace_all_in_children() failed"),
                $build_expr
            )
        }
    };
}

macro_rules! impl_binary_expr_common_methods {
    (
        $type_name:ident,
        $is_scalar:tt,
        $first_argument:ident,
        $second_argument:ident,
        $build_expr:expr,
        $with_substitute_zero_perturbations:tt
    ) => {
        impl_expr_common_methods!($is_scalar);

        #[inline]
        fn has_unperturbed_term(&self) -> bool {
            self.$first_argument.has_unperturbed_term()
                && self.$second_argument.has_unperturbed_term()
        }

        impl_binary_expr_common_methods!(
            @binary_expr_substitute_zero_perturbations
            $type_name,
            $is_scalar,
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
            impl_binary_expr_arg_operation!(
                self,
                $is_scalar,
                $first_argument,
                $second_argument,
                |first_arg: expr_arc_ref_ty!()| first_arg.eliminate(parameter.clone(), perturbations, min_order),
                |second_arg: expr_arc_ref_ty!()| second_arg.eliminate(parameter, perturbations, min_order),
                concat!(stringify!($type_name), "::eliminate() failed"),
                $build_expr
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
                return impl_binary_zero_expr!(self, $is_scalar);
            }

            impl_binary_expr_arg_operation!(
                self,
                $is_scalar,
                $first_argument,
                $second_argument,
                |arg: expr_arc_ref_ty!()| arg.remove_one(s),
                concat!(stringify!($type_name), "::remove_one() failed"),
                $build_expr
            )
        }

        #[inline]
        fn remove_all(&self, set: &expr_set_ty!()) -> expr_result_ty!() {
            if self.match_any_self(set, false) {
                return impl_binary_zero_expr!(self, $is_scalar);
            }

            impl_binary_expr_arg_operation!(
                self,
                $is_scalar,
                $first_argument,
                $second_argument,
                |arg: expr_arc_ref_ty!()| arg.remove_all(set),
                concat!(stringify!($type_name), "::remove_all() failed"),
                $build_expr
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

            let retained_first =
                self.$first_argument.retain_one(s, include_derivatives).map_err(|e| {
                    $crate::public::generic_expression_error(
                        concat!(
                            stringify!($type_name),
                            "::retain_one() failed for ",
                            stringify!($first_argument)
                        ),
                        self,
                        Some(::std::boxed::Box::new(e)),
                    )
                })?;

            let retained_second =
                self.$second_argument.retain_one(s, include_derivatives).map_err(|e| {
                    $crate::public::generic_expression_error(
                        concat!(
                            stringify!($type_name),
                            "::retain_one() failed for ",
                            stringify!($second_argument)
                        ),
                        self,
                        Some(::std::boxed::Box::new(e)),
                    )
                })?;

            let first_is_zero = $crate::public::is_zero_expr(&retained_first, None);
            let second_is_zero = $crate::public::is_zero_expr(&retained_second, None);

            if first_is_zero && second_is_zero {
                return impl_binary_zero_expr!(self, $is_scalar);
            }

            let new_first = if first_is_zero {
                self.$first_argument.clone()
            } else {
                retained_first
            };

            let new_second = if second_is_zero {
                self.$second_argument.clone()
            } else {
                retained_second
            };

            let first_changed = !::std::sync::Arc::ptr_eq(&new_first, &self.$first_argument)
                && &new_first != &self.$first_argument;

            let second_changed = !::std::sync::Arc::ptr_eq(&new_second, &self.$second_argument)
                && &new_second != &self.$second_argument;

            if !first_changed && !second_changed {
                Ok(self.clone_expr())
            } else {
                ($build_expr)(self, new_first, new_second)
            }
        }
    };

    (@binary_expr_substitute_zero_perturbations
        $type_name:ident,
        $is_scalar:tt,
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
            impl_binary_expr_arg_operation!(
                self,
                $is_scalar,
                $first_argument,
                $second_argument,
                |first_arg: expr_arc_ref_ty!()| first_arg.substitute_zero_perturbations(freq_tol.clone()),
                |second_arg: expr_arc_ref_ty!()| second_arg.substitute_zero_perturbations(freq_tol),
                concat!(stringify!($type_name), "::substitute_zero_perturbations() failed"),
                $build_expr
            )
        }
    };

    (@binary_expr_substitute_zero_perturbations
        $type_name:ident,
        $is_scalar:tt,
        $first_argument:ident,
        $second_argument:ident,
        $build_expr:expr,
        false
    ) => {};
}

macro_rules! impl_binary_zero_expr {
    ($self:ident, true) => {
        impl_zero_expr!(true)
    };

    ($self:ident, false) => {
        impl_zero_expr!(false)
    };

    // field case: called like `..., is_scalar, ...`
    ($self:ident, $field:ident) => {
        impl_zero_expr!($self.$field)
    };
}

macro_rules! impl_binary_expr_arg_operation {
    (
        $self:ident,
        $is_scalar:tt,
        $first_argument:ident,
        $second_argument:ident,
        $arg_operation:expr,
        $message:expr,
        $build_expr:expr
    ) => {{
        impl_binary_expr_arg_operation!(
            $self,
            $is_scalar,
            $first_argument,
            $second_argument,
            $arg_operation,
            $arg_operation,
            $message,
            $build_expr
        )
    }};

    (
        $self:ident,
        $is_scalar:tt,
        $first_argument:ident,
        $second_argument:ident,
        $first_arg_operation:expr,
        $second_arg_operation:expr,
        $message:expr,
        $build_expr:expr
    ) => {{
        let new_first = ($first_arg_operation)(&$self.$first_argument).map_err(|e| {
            $crate::public::generic_expression_error(
                concat!($message, " for ", stringify!($first_argument)),
                $self,
                Some(::std::boxed::Box::new(e)),
            )
        })?;

        if $crate::public::is_zero_expr(&new_first, None) {
            return impl_binary_zero_expr!($self, $is_scalar);
        }

        let new_second = ($second_arg_operation)(&$self.$second_argument).map_err(|e| {
            $crate::public::generic_expression_error(
                concat!($message, " for ", stringify!($second_argument)),
                $self,
                Some(::std::boxed::Box::new(e)),
            )
        })?;

        if $crate::public::is_zero_expr(&new_second, None) {
            return impl_binary_zero_expr!($self, $is_scalar);
        }

        if &new_first == &$self.$first_argument && &new_second == &$self.$second_argument {
            Ok($self.clone_expr())
        } else {
            ($build_expr)($self, new_first, new_second)
        }
    }};
}
