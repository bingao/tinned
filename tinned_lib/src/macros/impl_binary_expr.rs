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
        fn replace_expr_children(
            &self,
            map: &expr_map_ty!(),
            include_derivatives: bool,
        ) -> expr_result_ty!() {
            impl_binary_expr_arg_operation!(
                self,
                $is_scalar,
                $first_argument,
                $second_argument,
                |arg: expr_arc_ref_ty!()| arg.replace(map, include_derivatives),
                concat!(stringify!($type_name), "::replace_expr_children() failed"),
                $build_expr
            )
        }

        #[inline]
        fn retain_single(
            &self,
            s: expr_arc_ref_ty!(),
            include_derivatives: bool,
        ) -> expr_result_ty!() {
            if self.match_self_single(s, include_derivatives) {
                return Ok(self.clone_expr());
            }

            let retained_first =
                self.$first_argument.retain_single(s, include_derivatives).map_err(|e| {
                    $crate::public::generic_expression_error(
                        concat!(
                            stringify!($type_name),
                            "::retain_single() failed for ",
                            stringify!($first_argument)
                        ),
                        self,
                        Some(::std::boxed::Box::new(e)),
                    )
                })?;

            let retained_second =
                self.$second_argument.retain_single(s, include_derivatives).map_err(|e| {
                    $crate::public::generic_expression_error(
                        concat!(
                            stringify!($type_name),
                            "::retain_single() failed for ",
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
            parameter: expr_arc_ref_ty!(),
            perturbations: &[pert_arc_ty!()],
            min_order: u32,
        ) -> expr_result_ty!() {
            impl_binary_expr_arg_operation!(
                self,
                $is_scalar,
                $first_argument,
                $second_argument,
                |arg: expr_arc_ref_ty!()| arg.eliminate(parameter, perturbations, min_order),
                concat!(stringify!($type_name), "::eliminate() failed"),
                $build_expr
            )
        }

        #[inline]
        fn exist_any(&self, set: &expr_set_ty!(), include_derivatives: bool) -> bool {
            self.match_self_any(set, include_derivatives)
                || self.$first_argument.exist_any(set, include_derivatives)
                || self.$second_argument.exist_any(set, include_derivatives)
        }

        #[inline]
        fn find_superchains(&self, s: expr_arc_ref_ty!()) -> expr_differentiation_map_ty!() {
            if self.deep_eq_superchains(s) {
                ::std::collections::BTreeMap::from([(
                    self.total_order(),
                    ::std::collections::HashSet::from([self.clone_expr()]),
                )])
            } else {
                let mut result = self.$first_argument.find_superchains(s);

                for (order, subset) in self.$second_argument.find_superchains(s) {
                    result.entry(order).or_default().extend(subset);
                }

                result
            }
        }

        #[inline]
        fn remove(&self, set: &expr_set_ty!()) -> expr_result_ty!() {
            if self.match_self_any(set, false) {
                return impl_binary_zero_expr!(self, $is_scalar);
            }

            impl_binary_expr_arg_operation!(
                self,
                $is_scalar,
                $first_argument,
                $second_argument,
                |arg: expr_arc_ref_ty!()| arg.remove(set),
                concat!(stringify!($type_name), "::remove() failed"),
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
                |arg: expr_arc_ref_ty!()| arg.substitute_zero_perturbations(freq_tol.clone()),
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
        let new_first = ($arg_operation)(&$self.$first_argument).map_err(|e| {
            $crate::public::generic_expression_error(
                concat!($message, " for ", stringify!($first_argument)),
                $self,
                Some(::std::boxed::Box::new(e)),
            )
        })?;

        if $crate::public::is_zero_expr(&new_first, None) {
            return impl_binary_zero_expr!($self, $is_scalar);
        }

        let new_second = ($arg_operation)(&$self.$second_argument).map_err(|e| {
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
