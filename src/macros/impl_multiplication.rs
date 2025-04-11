macro_rules! impl_mul_traits {
    ($type_name:ident, $is_scalar:tt) => {
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
                    self == mul
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
                // Precompute the derivative of each factor and store it
                let diff_factors: Vec<Arc<dyn Expr>> =
                    self.factors.iter().map(|f| f.differentiate(s)).collect::<Result<_, _>>()?;

                let result = impl_mul_traits!(@finalize_differentiation self diff_factors s $is_scalar);

                result
            }
        }

        impl PartialEq for $type_name {
            fn eq(&self, other: &Self) -> bool {
                &self.coefficient == &other.coefficient && self.factors == other.factors
            }
        }

        impl Eq for $type_name {}

        impl std::fmt::Display for $type_name {
            // Format: coefficient * factor1 * factor2 * ..., omit coefficient if one
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                let wrote_coef = impl_mul_traits!(@non_one_coefficient self.coefficient, $is_scalar);

                if wrote_coef {
                    write!(f, "{}", self.coefficient)?;
                }

                if let Some((first, rest)) = self.factors().split_first() {
                    if wrote_coef {
                        write!(f, " * {}", first)?;
                    } else {
                        write!(f, "{}", first)?;
                    }

                    for factor in rest {
                        write!(f, " * {}", factor)?;
                    }
                }

                Ok(())
            }
        }
    };

    (@finalize_differentiation $self:ident $diff_factors:ident $s:ident true) => {{
        let mut results = Vec::new();

        for (i, diff) in $diff_factors.iter().enumerate() {
            // Skip derivative = 0 to avoid 0 * others = 0
            if is_zero_expr(diff) {
                continue;
            }

            let mut new_terms = $self.factors.clone();
            // For each factor i, replace it with its derivative while keeping
            // others intact
            new_terms[i] = diff.clone();
            new_terms.push($self.coefficient.clone().into());

            results.push(Self::new(new_terms)?);
        }

        crate::expressions::Add::new(results)
    }};

    (@finalize_differentiation $self:ident $diff_factors:ident $s:ident false) => {{
        let mut results = Vec::new();

        for (i, diff) in $diff_factors.iter().enumerate() {
            // Skip derivative = 0 to avoid 0 * others = 0
            if is_zero_expr(diff) {
                continue;
            }

            let mut new_terms = $self.factors.clone();
            // For each factor i, replace it with its derivative while keeping
            // others intact
            new_terms[i] = diff.clone();
            new_terms.push($self.coefficient.clone());

            results.push(Self::new(new_terms)?);
        }

        let diff_coef = $self.coefficient.differentiate($s)?;
        // If coefficient's derivative is non-zero, append it as one result
        if !is_zero_expr(&diff_coef) {
            let mut new_terms = $self.factors.clone();
            new_terms.push(diff_coef);
            results.push(Self::new(new_terms)?);
        }

        crate::expressions::MatrixAdd::new(results)
    }};

    (@non_one_coefficient $coefficient:expr, true) => { !$coefficient.is_one() };

    (@non_one_coefficient $coefficient:expr, false) => { !crate::utils::is_one_expr(&$coefficient) };
}
