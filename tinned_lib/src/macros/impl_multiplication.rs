macro_rules! impl_mul_traits {
    (
        $type_name:ident,
        $is_scalar:tt,
        $hash_delimiter:ident,
        $fmt_delimiter:ident
    ) => {
        impl $crate::core::ExprInternal for $type_name {
            impl_expr_internal_methods!($type_name, false);

            #[inline]
            fn hash_key(&self) -> String {
                format!(
                    "{}({}{}{})",
                    stringify!($type_name),
                    self.coefficient.hash_key(),
                    $hash_delimiter,
                    $crate::internal::join_mapped(&self.factors, $hash_delimiter, |f| f.hash_key())
                )
            }

            // For unambiguous replacement, we require equality for the
            // whole `Mul` so that we do not override methods
            // `eq_by_superchains()` and `replace_expr_self()` of
            // `ExprInternal`.

            fn replace_expr_fields(
                &self,
                map: &expr_map_ty!(),
                exact_equality: bool,
            ) -> expr_result_ty!() {
                impl_mul_traits!(
                    @mul_termwise_operation
                    self,
                    |factor: expr_arc_ref_ty!()| factor.replace(map, exact_equality),
                    concat!(stringify!($type_name), "::replace_expr_fields() failed"),
                    $is_scalar
                )
            }

            // FIXME: This method should be tested
            fn retain_expr_fields(
                &self,
                expr: &expr_arc_ty!(),
                exact_equality: bool,
            ) -> expr_result_ty!() {
                let mut num_changes = 0;
                // `bool` indicates if the new factor is different from the
                // original one or not
                let mut new_factors: ::std::vec::Vec<(expr_arc_ty!(), bool)>
                    = ::std::vec::Vec::with_capacity(self.factors.len());

                for factor in &self.factors {
                    let new_factor = factor.retain_expr(expr, exact_equality).map_err(|e| {
                        $crate::public::generic_expression_error(
                            concat!(stringify!($type_name), "::retain_expr_fields() failed"),
                            self,
                            Some(::std::boxed::Box::new(e)),
                        )
                    })?;
                    // `MatrixMul` or `Mul` will be retained as a whole if this
                    // factor retains completely.
                    if &new_factor == factor {
                        return Ok(self.clone_expr());
                    // This factor does not match any given ones. We save it in
                    // case that there is other factor(s) retained.
                    } else if $crate::public::is_zero_expr(&new_factor, None) {
                        new_factors.push((factor.clone(), false));
                    } else {
                        // Suppose `MatrixMul` is A*B*C*... = (Ak+Ar)*B*C*...,
                        // where Ak will be kept and Ar will be removed. The
                        // result after removal will be Ak*B*C*..., or
                        // A*B*C*... - Ar*B*C*..., where Ar = A - Ak.
                        new_factors.push((new_factor, true));
                        num_changes += 1;
                    }
                }

                let (new_coef, new_mul) = impl_mul_traits!(
                    @mul_coef_operation
                    self.coefficient,
                    |coef: expr_arc_ref_ty!()| coef.retain_expr(expr, exact_equality),
                    concat!(stringify!($type_name), "::retain_expr_fields() failed"),
                    $is_scalar
                );
                // Returns `MatrixMul` or `Mul` as a whole if the coefficient
                // retains completely.
                if !new_mul {
                    return Ok(self.clone_expr());
                }

                // The coefficient does not match any expression
                if $crate::public::is_zero_expr(&new_coef, None) {
                    match num_changes {
                        0 => impl_zero_expr!($is_scalar),
                        1 => {
                            // Only one factor retains partially, we simply
                            // return coefficient*Ak*B*C*...
                            let mut terms: expr_vec_ty!() = new_factors
                                .into_iter()
                                .map(|(factor, _changed)| factor)
                                .collect();
                            terms.push(
                                impl_mul_traits!(
                                    @mul_clone_coefficient
                                    self.coefficient,
                                    $is_scalar
                                )
                            );
                            Self::new(terms)
                        },
                        // As aforementioned, when there are factors partially
                        // retained, the result can be computed as
                        // A*B*C*...*R*S*T*... - Ar*Br*Cr*...*R*S*T*..., where
                        // Ar, Br, Cr, ... are parts that are removed, R, S, T,
                        // ... are those without retained parts.
                        _ => {
                            let mut terms: expr_vec_ty!() = new_factors
                                .into_iter()
                                .zip(self.factors.iter())
                                .map(|((factor, changed), original)| {
                                    if changed {
                                        $crate::public::subtract_exprs(original.clone(), factor)
                                    } else {
                                        Ok(factor)
                                    }
                                })
                                .collect::<::std::result::Result<::std::vec::Vec<_>, $crate::core::TinnedError>>()?;
                            terms.push(
                                impl_mul_traits!(
                                    @mul_clone_coefficient
                                    self.coefficient,
                                    $is_scalar
                                )
                            );
                            $crate::public::subtract_exprs(self.clone_expr(), Self::new(terms)?)
                        },
                    }
                // The coefficient retains partially
                } else {
                    match num_changes {
                        0 => {
                            // Only the coefficient retains partially, we
                            // return `new_coef`*A*B*C*...
                            let mut terms: expr_vec_ty!() = new_factors
                                .into_iter()
                                .map(|(factor, _changed)| factor)
                                .collect();
                            terms.push(new_coef);
                            Self::new(terms)
                        },
                        _ => {
                            let mut terms: expr_vec_ty!() = new_factors
                                .into_iter()
                                .zip(self.factors.iter())
                                .map(|((factor, changed), original)| {
                                    if changed {
                                        $crate::public::subtract_exprs(original.clone(), factor)
                                    } else {
                                        Ok(factor)
                                    }
                                })
                                .collect::<::std::result::Result<::std::vec::Vec<_>, $crate::core::TinnedError>>()?;
                            terms.push(
                                $crate::public::subtract_exprs(
                                    impl_mul_traits!(
                                        @mul_clone_coefficient
                                        self.coefficient,
                                        $is_scalar
                                    ),
                                    new_coef,
                                )?
                            );
                            $crate::public::subtract_exprs(self.clone_expr(), Self::new(terms)?)
                        },
                    }
                }
            }
        }

        #[::typetag::serde]
        impl $crate::core::Expr for $type_name {
            impl_expr_common_methods!($is_scalar);

            fn apply_zero_rules(
                &self,
                freq_tol: ::std::option::Option<$crate::public::NumberTolerance>,
            ) -> expr_result_ty!() {
                impl_mul_traits!(
                    @mul_termwise_operation
                    self,
                    |factor: expr_arc_ref_ty!()| factor.apply_zero_rules(freq_tol.clone()),
                    concat!(stringify!($type_name), "::apply_zero_rules() failed"),
                    $is_scalar
                )
            }

            fn differentiate(
                &self,
                s: &pert_arc_ty!(),
            ) -> expr_result_ty!() {
                // Precompute the derivative of each factor and store it
                let with_context = |f: expr_arc_ref_ty!()| {
                    f.differentiate(s).map_err(|e| {
                        $crate::public::generic_expression_error(
                            concat!(stringify!($type_name), "::differentiate() failed for factors"),
                            self,
                            Some(::std::boxed::Box::new(e)),
                        )
                    })
                };

                let diff_factors: expr_vec_ty!()
                    = self.factors.iter().map(with_context).collect::<::std::result::Result<_, _>>()?;

                impl_mul_traits!(
                    @mul_build_diff_expr
                    $type_name,
                    self,
                    diff_factors,
                    s,
                    $is_scalar
                )
            }

            fn eliminate(
                &self,
                parameter: &expr_arc_ty!(),
                perturbations: &[pert_arc_ty!()],
                min_order: u32,
            ) -> expr_result_ty!() {
                impl_mul_traits!(
                    @mul_termwise_operation
                    self,
                    |factor: expr_arc_ref_ty!()| {
                        factor.eliminate(parameter, perturbations, min_order)
                    },
                    concat!(stringify!($type_name), "::eliminate() failed"),
                    $is_scalar
                )
            }

            #[inline]
            fn exist_any(&self, set: &expr_set_ty!()) -> bool {
                if self.factors.iter().any(|factor| factor.exist_any(set)) {
                    return true;
                }

                set.iter().any(|expr| self.eq_expr(expr.as_ref()))
                    || self.coefficient.exist_any(set)
            }

            fn find_superchains(
                &self,
                s: &expr_arc_ty!(),
            ) -> expr_differentiation_map_ty!() {
                if self.deep_eq_superchains(s) {
                    return ::std::collections::BTreeMap::from([(
                        self.total_order(),
                        ::std::collections::HashSet::from([self.clone_expr()]),
                    )]);
                }

                let mut result: expr_differentiation_map_ty!() = ::std::collections::BTreeMap::new();
                for factor in &self.factors {
                    for (order, subset) in factor.find_superchains(s) {
                        result.entry(order).or_default().extend(subset);
                    }
                }

                if result.is_empty() {
                    return self.coefficient.find_superchains(s);
                }

                result
            }

            fn remove(&self, set: &expr_set_ty!()) -> expr_result_ty!() {
                if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                    return impl_zero_expr!($is_scalar);
                }

                impl_mul_traits!(
                    @mul_termwise_operation
                    self,
                    |factor: expr_arc_ref_ty!()| factor.remove(set),
                    concat!(stringify!($type_name), "::remove() failed"),
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

        impl ::std::fmt::Display for $type_name {
            // Format: coefficient * factor1 * factor2 * ..., omit coefficient if one
            fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
                let wrote_coef
                    = impl_mul_traits!(@mul_wrote_coefficient self.coefficient, $is_scalar);

                if wrote_coef {
                    write!(
                        f,
                        "{}{}{}",
                        self.coefficient,
                        $fmt_delimiter,
                        $crate::internal::join_mapped(&self.factors, $fmt_delimiter, |f| f.to_string()),
                    )
                } else {
                    write!(
                        f,
                        "{}",
                        $crate::internal::join_mapped(&self.factors, $fmt_delimiter, |f| f.to_string()),
                    )
                }
            }
        }
    };

    (@mul_build_diff_expr $type_name:ident, $self:ident, $diff_factors:ident, $s:ident, true) => {{
        let mut results = ::std::vec::Vec::with_capacity($diff_factors.len());

        for (i, diff) in $diff_factors.iter().enumerate() {
            // Skip derivative = 0 to avoid 0 * others = 0
            if $crate::public::is_zero_expr(diff, None) {
                continue;
            }

            //FIXME: collect common factors
            let mut new_terms = $self.factors.clone();
            // For each factor i, replace it with its derivative while keeping
            // others intact
            new_terms[i] = diff.clone();
            new_terms.push($self.coefficient.clone().into());

            results.push(Self::new(new_terms)?);
        }

        $crate::expressions::Add::new(results)
    }};

    (@mul_build_diff_expr $type_name:ident, $self:ident, $diff_factors:ident, $s:ident, false) => {{
        let mut results = ::std::vec::Vec::with_capacity($diff_factors.len() + 1);

        for (i, diff) in $diff_factors.iter().enumerate() {
            // Skip derivative = 0 to avoid 0 * others = 0
            if $crate::public::is_zero_expr(diff, None) {
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
                $crate::public::generic_expression_error(
                    concat!(stringify!($type_name), "::differentiate() failed for coefficient"),
                    $self,
                    Some(::std::boxed::Box::new(e)),
                )
            })?;
        // If coefficient's derivative is non-zero, append it as one result
        if !$crate::public::is_zero_expr(&diff_coef, None) {
            let mut new_terms = $self.factors.clone();
            new_terms.push(diff_coef);
            results.push(Self::new(new_terms)?);
        }

        $crate::expressions::MatrixAdd::new(results)
    }};

    (@mul_termwise_operation $self:ident, $operation:expr, $message:expr, $is_scalar:tt) => {{
        let (new_coef, mut new_mul) = impl_mul_traits!(
            @mul_coef_operation
            $self.coefficient,
            $operation,
            $message,
            $is_scalar
        );
        if $crate::public::is_zero_expr(&new_coef, None) {
            return impl_zero_expr!($is_scalar);
        }

        let mut new_factors = ::std::vec::Vec::with_capacity($self.factors.len() + 1);
        new_factors.push(new_coef);

        for factor in &$self.factors {
            let new_factor = ($operation)(factor).map_err(|e| {
                $crate::public::generic_expression_error(
                    concat!($message, " for factor"),
                    $self,
                    Some(::std::boxed::Box::new(e)),
                )
            })?;
            if $crate::public::is_zero_expr(&new_factor, None) {
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

    (@mul_wrote_coefficient $coefficient:expr, true) => {
        !$coefficient.is_one(None)
    };

    (@mul_wrote_coefficient $coefficient:expr, false) => {
        !$crate::public::is_one_expr(&$coefficient, None)
    };

    (@mul_coef_operation $coefficient:expr, $operation:expr, $message:expr, true) => {{
        let coef: expr_arc_ty!() = $coefficient.clone().into();
        let new_coef = ($operation)(&coef).map_err(|e| {
            $crate::public::generic_expression_error(
                concat!($message, " for coefficient"),
                &$coefficient,
                Some(::std::boxed::Box::new(e)),
            )
        })?;

        (new_coef.clone(), &new_coef != &coef)
    }};

    (@mul_coef_operation $coefficient:expr, $operation:expr, $message:expr, false) => {{
        let new_coef = ($operation)(&$coefficient).map_err(|e| {
            $crate::public::expression_error(
                concat!($message, " for coefficient"),
                &$coefficient,
                Some(::std::boxed::Box::new(e)),
            )
        })?;

        (new_coef.clone(), &new_coef != &$coefficient)
    }};

    (@mul_clone_coefficient $coefficient:expr, true) => {
        $coefficient.clone().into()
    };

    (@mul_clone_coefficient $coefficient:expr, false) => {
        $coefficient.clone()
    };
}
