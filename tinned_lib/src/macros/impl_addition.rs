macro_rules! impl_add_traits {
    ($type_name:ident, $hash_delimiter:ident, $fmt_delimiter:ident, $is_scalar:tt) => {
        impl ExprInternal for $type_name {
            impl_expr_internal_methods!($type_name, false);

            #[inline]
            fn hash_key(&self) -> String {
                format!(
                    "{}({})",
                    stringify!($type_name),
                    multi_expression_hash(&self.terms, $hash_delimiter),
                )
            }

            // For unambiguous replacement, we requirement equality for the
            // whole `Add` so that we do not override methods
            // `eq_by_superchains()` and `replace_expr_self()` of
            // `ExprInternal`.

            fn replace_expr_fields(
                &self,
                map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
                exact_equality: bool,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                impl_add_traits!(
                    @add_termwise_operation
                    self,
                    |term: &Arc<dyn Expr>| term.replace(map, exact_equality),
                    concat!(stringify!($type_name), "::replace_expr_fields() failed"),
                    $is_scalar
                )
            }

            fn retain_expr_fields(
                &self,
                set: &HashSet<Arc<dyn Expr>>,
                exact_equality: bool,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                impl_add_traits!(
                    @add_termwise_operation
                    self,
                    |term: &Arc<dyn Expr>| term.retain(set, exact_equality),
                    concat!(stringify!($type_name), "::retain_expr_fields() failed"),
                    $is_scalar
                )
            }
        }

        #[typetag::serde]
        impl Expr for $type_name {
            impl_expr_common_methods!($is_scalar);

            fn clean_temporum(
                &self,
                freq_tol: Option<NumberTolerance>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                impl_add_traits!(
                    @add_termwise_operation
                    self,
                    |term: &Arc<dyn Expr>| term.clean_temporum(freq_tol.clone()),
                    concat!(stringify!($type_name), "::clean_temporum() failed"),
                    $is_scalar
                )
            }

            fn differentiate(
                &self,
                s: &Arc<Perturbation>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                let mut diff_terms = Vec::with_capacity(self.terms.len());

                for term in &self.terms {
                    let diff = term.differentiate(s).map_err(|e| {
                        generic_expression_error(
                            concat!(stringify!($type_name), "::differentiate() failed"),
                            self,
                            Some(Box::new(e)),
                        )
                    })?;
                    if !is_zero_expr(&diff, None) {
                        diff_terms.push(diff);
                    }
                }

                Self::new(diff_terms)
            }

            fn eliminate(
                &self,
                parameter: &Arc<dyn Expr>,
                perturbations: &[Arc<Perturbation>],
                min_order: u32,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                impl_add_traits!(
                    @add_termwise_operation
                    self,
                    |term: &Arc<dyn Expr>| term.eliminate(parameter, perturbations, min_order),
                    concat!(stringify!($type_name), "::eliminate() failed"),
                    $is_scalar
                )
            }

            #[inline]
            fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
                if self.terms.iter().any(|term| term.exist_any(set)) {
                    return true;
                }

                set.iter().any(|expr| self.eq_expr(expr.as_ref()))
            }

            fn find_superchains(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
                if self.deep_eq_superchains(s) {
                    return BTreeMap::from([(self.total_order(), HashSet::from([self.clone_expr()]))]);
                }

                let mut result: BTreeMap<u32, HashSet<Arc<dyn Expr>>> = BTreeMap::new();
                for term in &self.terms {
                    for (order, subset) in term.find_superchains(s) {
                        result.entry(order).or_default().extend(subset);
                    }
                }

                result
            }

            fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
                if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                    return impl_zero_expr!($is_scalar);
                }

                impl_add_traits!(
                    @add_termwise_operation
                    self,
                    |term: &Arc<dyn Expr>| term.remove(set),
                    concat!(stringify!($type_name), "::remove() failed"),
                    $is_scalar
                )
            }
        }

        impl PartialEq for $type_name {
            fn eq(&self, other: &Self) -> bool {
                self.terms == other.terms
            }
        }

        impl Eq for $type_name {}

        impl std::fmt::Display for $type_name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "({})", multi_expression_format(&self.terms, $fmt_delimiter))
            }
        }
    };

    (@add_termwise_operation $self:ident, $operation:expr, $message:expr, $is_scalar:tt) => {{
        let mut new_terms = Vec::with_capacity($self.terms.len());
        let mut new_add = false;

        for term in &$self.terms {
            let new_term = ($operation)(term)
                .map_err(|e| generic_expression_error($message, $self, Some(Box::new(e))))?;
            if is_zero_expr(&new_term, None) {
                new_add = true;
            } else {
                if !new_add {
                    new_add = &new_term != term;
                }
                new_terms.push(new_term);
            }
        }

        if new_add {
            Self::new(new_terms)
        } else {
            Ok($self.clone_expr())
        }
    }};
}
