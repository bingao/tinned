macro_rules! impl_binary_expr_internal_methods {
    (
        $type_name:ident,
        $first_argument:ident,
        $second_argument:ident,
        $has_derivative:tt,
        $build_expr:expr
    ) => {
        impl_expr_internal_methods!($type_name, $has_derivative);

        #[inline]
        fn replace_expr_fields(
            &self,
            map: &expr_map_ty!(),
            exact_equality: bool,
        ) -> expr_result_ty!() {
            impl_binary_expr_arg_operation!(
                self,
                $first_argument,
                $second_argument,
                |arg: expr_arc_ref_ty!()| arg.replace(map, exact_equality),
                concat!(stringify!($type_name), "::replace_expr_fields() failed"),
                $build_expr
            )
        }

        #[inline]
        fn retain_expr_fields(
            &self,
            expr: expr_arc_ref_ty!(),
            exact_equality: bool,
        ) -> expr_result_ty!() {
            impl_binary_expr_arg_operation!(
                self,
                $first_argument,
                $second_argument,
                |arg: expr_arc_ref_ty!()| arg.retain_expr(expr, exact_equality),
                concat!(stringify!($type_name), "::retain_expr_fields() failed"),
                $build_expr
            )
        }
    };
}

macro_rules! impl_binary_expr_common_methods {
    (
        $type_name:ident,
        $first_argument:ident,
        $second_argument:ident,
        $is_scalar:tt,
        $build_expr:expr,
        $with_apply_zero_rules:tt
    ) => {
        impl_expr_common_methods!($is_scalar);

        impl_binary_expr_common_methods!(
            @binary_expr_apply_zero_rules
            $type_name,
            $first_argument,
            $second_argument,
            $build_expr,
            $with_apply_zero_rules
        );

        #[inline]
        fn eliminate(
            &self,
            parameter: expr_arc_ref_ty!(),
            perturbations: &[::std::sync::Arc<$crate::perturbations::Perturbation>],
            min_order: u32,
        ) -> expr_result_ty!() {
            impl_binary_expr_arg_operation!(
                self,
                $first_argument,
                $second_argument,
                |arg: expr_arc_ref_ty!()| arg.eliminate(parameter, perturbations, min_order),
                concat!(stringify!($type_name), "::eliminate() failed"),
                $build_expr
            )
        }

        #[inline]
        fn exist_any(&self, set: &expr_set_ty!()) -> bool {
            set.iter().any(|expr| self.eq_expr(expr.as_ref()))
                || self.$first_argument.exist_any(set)
                || self.$second_argument.exist_any(set)
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
            if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                impl_binary_expr_common_methods!(@binary_expr_return_zero self, $is_scalar);
            }

            impl_binary_expr_arg_operation!(
                self,
                $first_argument,
                $second_argument,
                |arg: expr_arc_ref_ty!()| arg.remove(set),
                concat!(stringify!($type_name), "::remove() failed"),
                $build_expr
            )
        }
    };

    (@binary_expr_apply_zero_rules
        $type_name:ident,
        $first_argument:ident,
        $second_argument:ident,
        $build_expr:expr,
        true
    ) => {
        #[inline]
        fn apply_zero_rules(
            &self,
            freq_tol: ::std::option::Option<$crate::public::NumberTolerance>,
        ) -> expr_result_ty!() {
            impl_binary_expr_arg_operation!(
                self,
                $first_argument,
                $second_argument,
                |arg: expr_arc_ref_ty!()| arg.apply_zero_rules(freq_tol.clone()),
                concat!(stringify!($type_name), "::apply_zero_rules() failed"),
                $build_expr
            )
        }
    };

    (@binary_expr_apply_zero_rules
        $type_name:ident,
        $first_argument:ident,
        $second_argument:ident,
        $build_expr:expr,
        false
    ) => {};

    (@binary_expr_return_zero $self:ident, true) => {
        return impl_zero_expr!(true);
    };

    (@binary_expr_return_zero $self:ident, false) => {
        return impl_zero_expr!(false);
    };

    // field case: called like `..., is_scalar, ...`
    (@binary_expr_return_zero $self:ident, $field:ident) => {
        return impl_zero_expr!($self.$field);
    };
}

macro_rules! impl_binary_expr_arg_operation {
    (
        $self:ident,
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

        let new_second = ($arg_operation)(&$self.$second_argument).map_err(|e| {
            $crate::public::generic_expression_error(
                concat!($message, " for ", stringify!($second_argument)),
                $self,
                Some(::std::boxed::Box::new(e)),
            )
        })?;

        if &new_first == &$self.$first_argument && &new_second == &$self.$second_argument {
            Ok($self.clone_expr())
        } else {
            ($build_expr)($self, new_first, new_second)
        }
    }};
}
