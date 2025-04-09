macro_rules! impl_mul_traits {
    ($type_name:ident, $is_scalar:literal) => {
        #[typetag::serde]
        impl Expr for $type_name {
            #[inline]
            fn as_any(&self) -> &dyn std::any::Any {
                self
            }

            #[inline]
            fn hash_key(&self) -> String {
                let keys: Vec<String> = self.factors.iter().map(|f| f.hash_key()).collect();
                format!(
                    "{}({}; {})",
                    stringify!($type_name),
                    self.coefficient.hash_key(),
                    keys.join(",")
                )
            }

            #[inline]
            fn is_scalar(&self) -> bool {
                $is_scalar
            }

            #[inline]
            fn eq_expr(&self, other: &dyn Expr) -> bool {
                if let Some(mul) = downcast_from_ref::<$type_name>(other) {
                    self.coefficient == mul.coefficient && self.factors == mul.factors
                } else {
                    false
                }
            }

            fn differentiate(
                &self,
                s: &crate::perturbations::Perturbation,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                // Precompute the derivative of each factor and store it
                let diff_factors: Vec<_> =
                    self.factors.iter().map(|f| f.differentiate(s)?).collect();

                let mut results = Vec::new();

                for (i, diff) in diff_factors.iter().enumerate() {
                    // Skip derivative = 0 to avoid 0 * others = 0
                    if is_zero_expr(diff) {
                        continue;
                    }

                    let mut new_terms = self.factors.clone();
                    // For each factor i, replace it with its derivative while keeping
                    // others intact
                    new_terms[i] = diff.clone();
                    if $is_scalar {
                        new_terms.push(self.coefficient.into());
                    } else {
                        new_terms.push(self.coefficient.clone());
                    }

                    results.push(Self::new(new_terms)?);
                }

                if $is_scalar {
                    crate::expressions::Add::new(results)
                } else {
                    let diff_coef = self.coefficient.differentiate(s)?;
                    // If coefficient's derivative is non-zero, append it as one result
                    if !is_zero_expr(&diff_coef) {
                        let mut new_terms = self.factors.clone();
                        new_terms.push(diff_coef);
                        results.push(Self::new(new_terms)?);
                    }

                    crate::expressions::MatrixAdd::new(results)
                }
            }
        }

        impl std::fmt::Display for $type_name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                if $is_scalar {
                    crate::utils::fmt_mul(f, self)
                } else {
                    crate::utils::fmt_matrix_mul(f, self)
                }
            }
        }
    };
}
