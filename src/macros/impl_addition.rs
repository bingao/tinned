macro_rules! impl_add_traits {
    ($type_name:ident, $is_scalar:literal) => {
        #[typetag::serde]
        impl Expr for $type_name {
            #[inline]
            fn as_any(&self) -> &dyn std::any::Any {
                self
            }

            #[inline]
            fn hash_key(&self) -> String {
                let keys: Vec<String> = self.terms.iter().map(|t| t.hash_key()).collect();
                format!("{}({})", stringify!($type_name), keys.join("+"))
            }

            #[inline]
            fn is_scalar(&self) -> bool {
                $is_scalar
            }

            #[inline]
            fn eq_expr(&self, other: &dyn Expr) -> bool {
                if let Some(add) = downcast_from_ref::<$type_name>(other) {
                    self.terms == add.terms
                } else {
                    false
                }
            }

            fn differentiate(
                &self,
                s: &crate::perturbations::Perturbation,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                let mut diff_terms = Vec::new();

                for term in &self.terms {
                    let diff = term.differentiate(s)?;
                    if !crate::utils::is_zero_expr(&diff) {
                        diff_terms.push(diff);
                    }
                }

                Self::new(diff_terms)
            }
        }

        impl std::fmt::Display for $type_name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "(")?;
                let mut iter = self.terms.iter();
                if let Some(first) = iter.next() {
                    write!(f, "{}", first)?;
                    for term in iter {
                        write!(f, " + {}", term)?;
                    }
                }
                write!(f, ")")
            }
        }
    };
}
