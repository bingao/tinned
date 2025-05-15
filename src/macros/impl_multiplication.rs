macro_rules! impl_mul_traits {
    (
        $type_name:ident,
        $is_scalar:tt,
        $hash_delimiter:ident,
        $fmt_delimiter:ident
    ) => {
        impl ExprInternal for $type_name {
            impl_expr_internal_methods!($type_name);

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
        }

        #[typetag::serde]
        impl Expr for $type_name {
            impl_expr_common_methods!($is_scalar);

            fn clean_temporum(
                &self,
                freq_tol: Option<NumberTolerance>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                impl_mul_traits!(
                    @mul_termwise_operation
                    self,
                    |factor: &Arc<dyn Expr>| factor.clean_temporum(freq_tol.clone()),
                    concat!(stringify!($type_name), "::clean_temporum() failed"),
                    $is_scalar
                )
            }

            fn differentiate(
                &self,
                s: &Arc<Perturbation>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                // Precompute the derivative of each factor and store it
                let with_context = |f: &Arc<dyn Expr>| {
                    f.differentiate(s).map_err(|e| {
                        generic_expression_error(
                            concat!(stringify!($type_name), "::differentiate() failed for factors"),
                            self,
                            Some(Box::new(e)),
                        )
                    })
                };

                let diff_factors: Vec<Arc<dyn Expr>> =
                    self.factors.iter().map(with_context).collect::<Result<_, _>>()?;

                let result
                    = impl_mul_traits!(@build_differentiation_expr self, diff_factors, s, $is_scalar);

                result
            }

            fn eliminate(
                &self,
                parameter: &Arc<dyn Expr>,
                perturbations: &[Arc<Perturbation>],
                min_order: u32,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                impl_mul_traits!(
                    @mul_termwise_operation
                    self,
                    |factor: &Arc<dyn Expr>| factor.eliminate(parameter, perturbations, min_order),
                    concat!(stringify!($type_name), "::eliminate() failed"),
                    $is_scalar
                )
            }

            #[inline]
            fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
                if self.factors.iter().any(|factor| factor.exist_any(set)) {
                    return true;
                }

                set.iter().any(|expr| self.eq_expr(expr.as_ref()))
                    || self.coefficient.exist_any(set)
            }

            fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
                if self.match_for_find_all(s) {
                    return BTreeMap::from([(self.total_order(), HashSet::from([self.clone_expr()]))]);
                }

                let mut result: BTreeMap<u32, HashSet<Arc<dyn Expr>>> = BTreeMap::new();
                for factor in &self.factors {
                    for (order, subset) in factor.find_all(s) {
                        result.entry(order).or_default().extend(subset);
                    }
                }

               if result.is_empty() {
                   return self.coefficient.find_all(s);
               }

               result
            }

            fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
                if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                    return impl_zero_expr!($is_scalar);
                }

                impl_mul_traits!(
                    @mul_termwise_operation
                    self,
                    |factor: &Arc<dyn Expr>| factor.remove(set),
                    concat!(stringify!($type_name), "::remove() failed"),
                    $is_scalar
                )
            }

            fn replace(
                &self,
                map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                if let Some((_, value)) = map.iter().find(|(key, _)| self.eq_expr(key.as_ref())) {
                    return Ok(value.clone());
                }

                impl_mul_traits!(
                    @mul_termwise_operation
                    self,
                    |factor: &Arc<dyn Expr>| factor.replace(map),
                    concat!(stringify!($type_name), "::replace() failed"),
                    $is_scalar
                )
            }

            fn replace_all(
                &self,
                map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                // For unambiguous replacement, we requirement equality for the
                // whole `Mul`.
                if let Some((_, value)) = map.iter().find(|(key, _)| self.match_for_replace_all(key)) {
                    return Ok(value.clone());
                }

                impl_mul_traits!(
                    @mul_termwise_operation
                    self,
                    |factor: &Arc<dyn Expr>| factor.replace_all(map),
                    concat!(stringify!($type_name), "::replace_all() failed"),
                    $is_scalar
                )
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

    (@build_differentiation_expr $self:ident, $diff_factors:ident, $s:ident, true) => {{
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

    (@build_differentiation_expr $self:ident, $diff_factors:ident, $s:ident, false) => {{
        let mut results = Vec::with_capacity($diff_factors.len() + 1);

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
                    concat!(stringify!($type_name), "::differentiate() failed for coefficient"),
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

    (@mul_termwise_operation $self:ident, $operation:expr, $message:expr, $is_scalar:tt) => {{
        let mut new_factors = Vec::with_capacity($self.factors.len() + 1);
        impl_mul_traits!(@push_coefficient new_factors, $self.coefficient, $is_scalar);
        let mut new_mul = false;

        for factor in &$self.factors {
            let new_factor = ($operation)(factor)
                .map_err(|e| generic_expression_error($message, $self, Some(Box::new(e))))?;
            if is_zero_expr(&new_factor, None) {
                return impl_zero_expr!($is_scalar);
            } else {
                if !new_mul {
                    new_mul = &new_factor != factor;
                }
                new_factors.push(new_factor);
            }
        }

        if new_mul {
            Self::new(new_factors)
        } else {
            Ok($self.clone_expr())
        }
    }};

    (@non_one_coefficient $coefficient:expr, true) => { !$coefficient.is_one(None) };

    (@non_one_coefficient $coefficient:expr, false) => { !is_one_expr(&$coefficient, None) };

    (@push_coefficient $factors:expr, $coefficient:expr, true) => {
        $factors.push($coefficient.clone().into())
    };

    (@push_coefficient $factors:expr, $coefficient:expr, false) => {
        $factors.push($coefficient.clone())
    };
}
