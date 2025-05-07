macro_rules! impl_mul_traits {
    (
        $type_name:ident,
        $is_scalar:tt,
        $build_zero_expr:ident,
        $hash_delimiter:ident,
        $fmt_delimiter:ident
    ) => {
        #[typetag::serde]
        impl Expr for $type_name {
            #[inline]
            fn as_any(&self) -> &dyn std::any::Any {
                self
            }

            #[inline]
            fn hash_key(&self) -> String {
                format!(
                    "{}({}{}{})",
                    stringify!($type_name),
                    self.coefficient.hash_key(),
                    $hash_delimiter,
                    multi_expression_hash(&self.factors, $hash_delimiter),
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

            fn clean_temporum(
                &self,
                num_tol: Option<NumberTolerance>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                let mut new_factors = Vec::with_capacity(self.factors.len() + 1);
                new_factors.push(
                    impl_mul_traits!(@clone_coefficient self.coefficient, $is_scalar)
                );
                let mut new_mul = false;

                for factor in &self.factors {
                    let new_factor = factor.clean_temporum(num_tol).map_err(|e| {
                        generic_expression_error("clean_temporum() failed", self, Some(Box::new(e)))
                    })?;
                    if is_zero_expr(&new_factor, num_tol) {
                        return Ok($build_zero_expr());
                    } else {
                        if !new_mul {
                            new_mul = new_factor != factor;
                        }
                        new_factors.push(new_factor);
                    }
                }

                if new_mul {
                    Self::new(new_factors)
                } else {
                    Ok(Arc::new(self.clone_expr()))
                }
            }

            fn differentiate(
                &self,
                s: &Arc<crate::perturbations::Perturbation>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                // Precompute the derivative of each factor and store it
                let with_context = |f: &Arc<dyn Expr>| {
                    f.differentiate(s).map_err(|e| {
                        generic_expression_error("differentiate() failed", self, Some(Box::new(e)))
                    })
                };

                let diff_factors: Vec<Arc<dyn Expr>> =
                    self.factors.iter().map(with_context).collect::<Result<_, _>>()?;

                let result
                    = impl_mul_traits!(@finalize_differentiation self diff_factors s $is_scalar);

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
                let wrote_coef
                    = impl_mul_traits!(@non_one_coefficient self.coefficient, $is_scalar);

                if wrote_coef {
                    write!(
                        f,
                        "{}{}{}",
                        self.coefficient,
                        $fmt_delimiter,
                        multi_expression_format(&self.factors, $fmt_delimiter),
                    )
                } else {
                    write!(f, "{}", multi_expression_format(&self.factors, $fmt_delimiter))
                }
            }
        }
    };

    (@finalize_differentiation $self:ident $diff_factors:ident $s:ident true) => {{
        let mut results = Vec::with_capacity($diff_factors.len());

        for (i, diff) in $diff_factors.iter().enumerate() {
            // Skip derivative = 0 to avoid 0 * others = 0
            if is_zero_expr(diff, None) {
                continue;
            }

            let mut new_terms = $self.factors.clone();
            // For each factor i, replace it with its derivative while keeping
            // others intact
            new_terms[i] = diff.clone();
            new_terms.push($self.coefficient.clone().into());

            results.push(Self::new(new_terms)?);
        }

        Add::new(results)
    }};

    (@finalize_differentiation $self:ident $diff_factors:ident $s:ident false) => {{
        let mut results = Vec::with_capacity($diff_factors.len());

        for (i, diff) in $diff_factors.iter().enumerate() {
            // Skip derivative = 0 to avoid 0 * others = 0
            if is_zero_expr(diff, None) {
                continue;
            }

            let mut new_terms = $self.factors.clone();
            // For each factor i, replace it with its derivative while keeping
            // others intact
            new_terms[i] = diff.clone();
            new_terms.push($self.coefficient.clone());

            results.push(Self::new(new_terms)?);
        }

        let diff_coef = $self
            .coefficient
            .differentiate($s)
            .map_err(|e| {
                generic_expression_error(
                    "differentiate() on coefficient failed",
                    $self,
                    Some(Box::new(e)),
                )
            })?;
        // If coefficient's derivative is non-zero, append it as one result
        if !is_zero_expr(&diff_coef, None) {
            let mut new_terms = $self.factors.clone();
            new_terms.push(diff_coef);
            results.push(Self::new(new_terms)?);
        }

        MatrixAdd::new(results)
    }};

    (@non_one_coefficient $coefficient:expr, true) => { !$coefficient.is_one(None) };

    (@non_one_coefficient $coefficient:expr, false) => { !is_one_expr(&$coefficient, None) };

    (@clone_coefficient $coefficient:expr, true) => { $coefficient.into() };

    (@clone_coefficient $coefficient:expr, false) => { $coefficient.clone() };
}
