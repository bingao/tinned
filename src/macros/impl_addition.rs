macro_rules! impl_add_traits {
    ($type_name:ident, $hash_delimiter:ident, $fmt_delimiter:ident, $is_scalar:literal) => {
        #[typetag::serde]
        impl Expr for $type_name {
            #[inline]
            fn as_any(&self) -> &dyn std::any::Any {
                self
            }

            #[inline]
            fn hash_key(&self) -> String {
                format!(
                    "{}({})",
                    stringify!($type_name),
                    multi_expression_hash(&self.terms, $hash_delimiter),
                )
            }

            #[inline]
            fn is_scalar(&self) -> bool {
                $is_scalar
            }

            #[inline]
            fn clone_expr(&self) -> Arc<dyn Expr> {
                Arc::new(self.clone())
            }

            #[inline]
            fn eq_expr(&self, other: &dyn Expr) -> bool {
                if let Some(add) = downcast_from_ref::<$type_name>(other) {
                    self == add
                } else {
                    false
                }
            }

            #[inline]
            fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "{self}")
            }

            fn clean_temporum(
                &self,
                freq_tol: Option<NumberTolerance>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                impl_add_traits!(
                    @add_termwise_operation
                    self,
                    |term: &Arc<dyn Expr>| term.clean_temporum(freq_tol.clone()),
                    concat!(stringify!($type_name), "clean_temporum() failed")
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
                            concat!(stringify!($type_name), "differentiate() failed"),
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
                    concat!(stringify!($type_name), "eliminate() failed")
                )
            }

            #[inline]
            fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
                if self.terms.iter().any(|term| term.exist_any(set)) {
                    return true;
                }

                set.iter().any(|expr| self.eq_expr(expr.as_ref()))
            }

            fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
                if self.eq_shallow(s) {
                    return BTreeMap::from([(self.total_order(), HashSet::from([self.clone_expr()]))]);
                }

                let mut result: BTreeMap<u32, HashSet<Arc<dyn Expr>>> = BTreeMap::new();
                for term in &self.terms {
                    for (order, subset) in term.find_all(s) {
                        result.entry(order).or_default().extend(subset);
                    }
                }

                result
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

    (@add_termwise_operation $self:ident, $operation:expr, $message:expr) => {{
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
