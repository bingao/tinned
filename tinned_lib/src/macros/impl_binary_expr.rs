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
            map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
            exact_equality: bool,
        ) -> Result<Arc<dyn Expr>, TinnedError> {
            impl_binary_expr_arg_operation!(
                self,
                $first_argument,
                $second_argument,
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
            impl_binary_expr_arg_operation!(
                self,
                $first_argument,
                $second_argument,
                |arg: &Arc<dyn Expr>| arg.retain(set, exact_equality),
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
        $with_clean_temporum:tt
    ) => {
        impl_expr_common_methods!($is_scalar);

        impl_binary_expr_common_methods!(
            @binary_expr_clean_temporum
            $type_name,
            $first_argument,
            $second_argument,
            $build_expr,
            $with_clean_temporum
        );

        #[inline]
        fn eliminate(
            &self,
            parameter: &Arc<dyn Expr>,
            perturbations: &[Arc<Perturbation>],
            min_order: u32,
        ) -> Result<Arc<dyn Expr>, TinnedError> {
            impl_binary_expr_arg_operation!(
                self,
                $first_argument,
                $second_argument,
                |arg: &Arc<dyn Expr>| arg.eliminate(parameter, perturbations, min_order),
                concat!(stringify!($type_name), "::eliminate() failed"),
                $build_expr
            )
        }

        #[inline]
        fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
            set.iter().any(|expr| self.eq_expr(expr.as_ref()))
                || self.$first_argument.exist_any(set)
                || self.$second_argument.exist_any(set)
        }

        #[inline]
        fn find_superchains(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
            if self.deep_eq_superchains(s) {
                BTreeMap::from([(self.total_order(), HashSet::from([self.clone_expr()]))])
            } else {
                let mut result = self.$first_argument.find_superchains(s);
                for (order, subset) in self.$second_argument.find_superchains(s) {
                    result.entry(order).or_default().extend(subset);
                }

                result
            }
        }

        #[inline]
        fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
            if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                impl_binary_expr_common_methods!(@binary_expr_return_zero self, $is_scalar);
            }

            impl_binary_expr_arg_operation!(
                self,
                $first_argument,
                $second_argument,
                |arg: &Arc<dyn Expr>| arg.remove(set),
                concat!(stringify!($type_name), "::remove() failed"),
                $build_expr
            )
        }
    };

    (@binary_expr_clean_temporum
        $type_name:ident,
        $first_argument:ident,
        $second_argument:ident,
        $build_expr:expr,
        true
    ) => {
        #[inline]
        fn clean_temporum(
            &self,
            freq_tol: Option<NumberTolerance>,
        ) -> Result<Arc<dyn Expr>, TinnedError> {
            impl_binary_expr_arg_operation!(
                self,
                $first_argument,
                $second_argument,
                |arg: &Arc<dyn Expr>| arg.clean_temporum(freq_tol.clone()),
                concat!(stringify!($type_name), "::clean_temporum() failed"),
                $build_expr
            )
        }
    };

    (@binary_expr_clean_temporum
        $type_name:ident,
        $first_argument:ident,
        $second_argument:ident,
        $build_expr:expr,
        false
    ) => { };

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
            generic_expression_error(
                concat!($message, " for ", stringify!($first_argument)),
                $self,
                Some(Box::new(e)),
            )
        })?;
        let new_second = ($arg_operation)(&$self.$second_argument).map_err(|e| {
            generic_expression_error(
                concat!($message, " for ", stringify!($second_argument)),
                $self,
                Some(Box::new(e)),
            )
        })?;

        if &new_first == &$self.$first_argument && &new_second == &$self.$second_argument {
            Ok($self.clone_expr())
        } else {
            ($build_expr)($self, new_first, new_second)
        }
    }};
}
