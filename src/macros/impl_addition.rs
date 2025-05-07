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
            fn clone_expr(&self) -> Self {
                self.clone()
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
                num_tol: Option<NumberTolerance>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                let mut new_terms = Vec::with_capacity(self.terms.len());
                let mut new_add = false;

                for term in &self.terms {
                    let new_term = term.clean_temporum(num_tol).map_err(|e| {
                        generic_expression_error("clean_temporum() failed", self, Some(Box::new(e)))
                    })?;
                    if is_zero_expr(&new_term, num_tol) {
                        new_add = true;
                    } else {
                        if !new_add {
                            new_add = new_term != term;
                        }
                        new_terms.push(new_term);
                    }
                }

                if new_add {
                    Self::new(new_terms)
                } else {
                    Ok(Arc::new(self.clone_expr()))
                }
            }

            fn differentiate(
                &self,
                s: &Arc<crate::perturbations::Perturbation>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                let mut diff_terms = Vec::with_capacity(self.terms.len());

                for term in &self.terms {
                    let diff = term.differentiate(s).map_err(|e| {
                        generic_expression_error("differentiate() failed", self, Some(Box::new(e)))
                    })?;
                    if !is_zero_expr(&diff, None) {
                        diff_terms.push(diff);
                    }
                }

                Self::new(diff_terms)
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
}
