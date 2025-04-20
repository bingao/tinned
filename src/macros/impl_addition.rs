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
                    join_exprs_for_hash(&self.terms, $hash_delimiter),
                )
            }

            #[inline]
            fn is_scalar(&self) -> bool {
                $is_scalar
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

            fn differentiate(
                &self,
                s: &Arc<crate::perturbations::Perturbation>,
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

        impl PartialEq for $type_name {
            fn eq(&self, other: &Self) -> bool {
                self.terms == other.terms
            }
        }

        impl Eq for $type_name {}

        impl std::fmt::Display for $type_name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "({})", join_exprs_for_display(&self.terms, $fmt_delimiter))
            }
        }
    };
}
