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
                $build_expr,
            )
        }

        #[inline]
        fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
            set.iter().any(|expr| self.eq_expr(expr.as_ref()))
                || self.$first_argument.exist_any(set)
                || self.$second_argument.exist_any(set)
        }

        #[inline]
        fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
            if self.match_for_find_all(s) {
                BTreeMap::from([(self.find_all_key(), HashSet::from([self.clone_expr()]))])
            } else {
                let mut result = self.$first_argument.find_all(s);
                for (order, subset) in self.$second_argument.find_all(s) {
                    result.entry(order).or_default().extend(subset);
                }

                result
            }
        }

        #[inline]
        fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
            if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                return impl_zero_expr!($is_scalar);
            }

            impl_binary_expr_arg_operation!(
                self,
                $first_argument,
                $second_argument,
                |arg: &Arc<dyn Expr>| arg.remove(set),
                concat!(stringify!($type_name), "::remove() failed"),
                $build_expr,
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

            impl_binary_expr_arg_operation!(
                self,
                $first_argument,
                $second_argument,
                |arg: &Arc<dyn Expr>| arg.replace(map),
                concat!(stringify!($type_name), "::replace() failed"),
                $build_expr,
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
                $build_expr,
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
}

macro_rules! impl_binary_expr_arg_operation {
    (
        $self:ident,
        $first_argument:ident,
        $second_argument:ident,
        $arg_operation:expr,
        $message:expr,
        $build_expr:expr,
    ) => {{
        let new_first = ($arg_operation)(&$self.$first_argument).map_err(|e| {
            generic_expression_error(
                concat!($message, " for first argument"),
                $self,
                Some(Box::new(e)),
            )
        })?;
        let new_second = ($arg_operation)(&$self.$second_argument).map_err(|e| {
            generic_expression_error(
                concat!($message, " for second argument"),
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
