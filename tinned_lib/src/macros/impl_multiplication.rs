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

            #[inline]
            fn expr_order(&self) -> u32 {
                let mut expr_order = self.coefficient.expr_order();

                for factor in &self.factors {
                    expr_order += factor.expr_order()
                }

                expr_order
            }

            // For unambiguous replacement, we require equality for the
            // whole `Mul` so that we do not override methods
            // `eq_by_superchains()` and `apply_replacement()` of
            // `ExprInternal`.

            fn replace_one_in_children(
                &self,
                expr: expr_arc_ref_ty!(),
                replacement: expr_arc_ty!(),
                include_derivatives: bool,
            ) -> expr_result_ty!() {
                impl_mul_traits!(
                    @mul_termwise_operation_any_zero
                    self,
                    |factor: expr_arc_ref_ty!()| factor.replace_one(expr, replacement.clone(), include_derivatives),
                    concat!(stringify!($type_name), "::replace_one_in_children() failed"),
                    $is_scalar
                )
            }

            fn replace_all_in_children(
                &self,
                map: &expr_map_ty!(),
                include_derivatives: bool,
            ) -> expr_result_ty!() {
                impl_mul_traits!(
                    @mul_termwise_operation_any_zero
                    self,
                    |factor: expr_arc_ref_ty!()| factor.replace_all(map, include_derivatives),
                    concat!(stringify!($type_name), "::replace_all_in_children() failed"),
                    $is_scalar
                )
            }
        }

        #[::typetag::serde]
        impl $crate::core::Expr for $type_name {
            impl_expr_common_methods!($is_scalar);

            #[inline]
            fn has_unperturbed_term(&self) -> bool {
                self.coefficient.has_unperturbed_term()
                    && self.factors.iter().all(|factor| factor.has_unperturbed_term())
            }

            fn substitute_zero_perturbations(
                &self,
                freq_tol: ::std::option::Option<$crate::public::NumberTolerance>,
            ) -> expr_result_ty!() {
                impl_mul_traits!(
                    @mul_termwise_operation_any_zero
                    self,
                    |factor: expr_arc_ref_ty!()| factor.substitute_zero_perturbations(freq_tol.clone()),
                    concat!(stringify!($type_name), "::substitute_zero_perturbations() failed"),
                    $is_scalar
                )
            }

            impl_mul_traits!(
                @impl_mul_differentiate
                $type_name,
                $is_scalar
            );

            fn eliminate(
                &self,
                parameter: expr_arc_ty!(),
                perturbations: &[pert_arc_ty!()],
                min_order: u32,
            ) -> expr_result_ty!() {
                impl_mul_traits!(
                    @mul_termwise_operation_any_zero
                    self,
                    |factor: expr_arc_ref_ty!()| {
                        factor.eliminate(parameter.clone(), perturbations, min_order)
                    },
                    concat!(stringify!($type_name), "::eliminate() failed"),
                    $is_scalar
                )
            }

            fn find_all(
                &self,
                s: expr_arc_ref_ty!(),
            ) -> expr_differentiation_map_ty!() {
                if self.deep_eq_superchains(s) {
                    return ::std::collections::BTreeMap::from([(
                        self.expr_order(),
                        ::std::collections::HashSet::from([self.clone_expr()]),
                    )]);
                }

                let mut result: expr_differentiation_map_ty!() = ::std::collections::BTreeMap::new();
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

            #[inline]
            fn match_one(&self, s: expr_arc_ref_ty!(), include_derivatives: bool) -> bool {
                self.match_one_self(s, include_derivatives)
                    || self.coefficient.match_one(s, include_derivatives)
                    || self.factors.iter().any(|factor| factor.match_one(s, include_derivatives))
            }

            #[inline]
            fn match_any(&self, set: &expr_set_ty!(), include_derivatives: bool) -> bool {
                self.match_any_self(set, include_derivatives)
                    || self.coefficient.match_any(set, include_derivatives)
                    || self.factors.iter().any(|factor| factor.match_any(set, include_derivatives))
            }

            fn remove_one(&self, s: expr_arc_ref_ty!()) -> expr_result_ty!() {
                if self.match_one_self(s, false) {
                    return impl_zero_expr!($is_scalar);
                }

                impl_mul_traits!(
                    @mul_termwise_operation_any_zero
                    self,
                    |factor: expr_arc_ref_ty!()| factor.remove_one(s),
                    concat!(stringify!($type_name), "::remove_one() failed"),
                    $is_scalar
                )
            }

            fn remove_all(&self, set: &expr_set_ty!()) -> expr_result_ty!() {
                if self.match_any_self(set, false) {
                    return impl_zero_expr!($is_scalar);
                }

                impl_mul_traits!(
                    @mul_termwise_operation_any_zero
                    self,
                    |factor: expr_arc_ref_ty!()| factor.remove_all(set),
                    concat!(stringify!($type_name), "::remove_all() failed"),
                    $is_scalar
                )
            }

            fn retain_one(
                &self,
                s: expr_arc_ref_ty!(),
                include_derivatives: bool,
            ) -> expr_result_ty!() {
                if self.match_one_self(s, include_derivatives) {
                    return Ok(self.clone_expr());
                }

                impl_mul_traits!(
                    @mul_termwise_operation_all_zero
                    self,
                    |factor: expr_arc_ref_ty!()| factor.retain_one(s, include_derivatives),
                    concat!(stringify!($type_name), "::retain_one() failed"),
                    $is_scalar
                )
            }

            fn retain_any(
                &self,
                set: &expr_set_ty!(),
                include_derivatives: bool,
            ) -> expr_result_ty!() {
                if self.match_any_self(set, include_derivatives) {
                    return Ok(self.clone_expr());
                }

                impl_mul_traits!(
                    @mul_termwise_operation_all_zero
                    self,
                    |factor: expr_arc_ref_ty!()| factor.retain_any(set, include_derivatives),
                    concat!(stringify!($type_name), "::retain_any() failed"),
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

    (@impl_mul_differentiate $type_name:ident, true) => {
        #[inline]
        fn differentiate(&self, s: pert_arc_ty!()) -> expr_result_ty!() {
            //FIXME: how to collect common factors?
            let results = $crate::internal::differentiate_operands(
                &self.factors,
                |term| term.differentiate(s.clone()),
                |mut new_factors| {
                    new_factors.push(self.coefficient.clone().into());
                    Self::new(new_factors)
                },
            )
            .map_err(|e| {
                $crate::public::generic_expression_error(
                    "Mul::differentiate() failed",
                    self,
                    Some(::std::boxed::Box::new(e)),
                )
            })?;

            $crate::expressions::Add::new(results)
        }
    };

    (@impl_mul_differentiate $type_name:ident, false) => {
        #[inline]
        fn differentiate(&self, s: pert_arc_ty!()) -> expr_result_ty!() {
            let results = $crate::internal::differentiate_operands_and_base(
                &self.factors,
                &self.coefficient,
                |term| term.differentiate(s.clone()),
                //FIXME: how to collect common factors?
                |mut new_factors| {
                    new_factors.push(self.coefficient.clone());
                    Self::new(new_factors)
                },
                |diff_coef| {
                    let mut new_factors = self.factors.clone();
                    new_factors.push(diff_coef);
                    Self::new(new_factors)
                },
            )
            .map_err(|e| {
                $crate::public::generic_expression_error(
                    "MatrixMul::differentiate() failed",
                    self,
                    Some(::std::boxed::Box::new(e)),
                )
            })?;

            $crate::expressions::MatrixAdd::new(results)
        }
    };

    // Termwise operation for `Mul`
    (@mul_termwise_operation_any_zero $self:ident, $operation:expr, $message:expr, true) => {{
        $crate::internal::transform_mul_any_zero(
            $self,
            &$self.coefficient.clone().into(),
            &$self.factors,
            |factor: expr_arc_ref_ty!()| ($operation)(factor),
            $message,
            |new_factors| Self::new(new_factors),
            || Ok($crate::expressions::Number::zero()),
        )
    }};

    // Termwise operation for `MatrixMul`
    (@mul_termwise_operation_any_zero $self:ident, $operation:expr, $message:expr, false) => {{
        $crate::internal::transform_mul_any_zero(
            $self,
            &$self.coefficient,
            &$self.factors,
            |factor: expr_arc_ref_ty!()| ($operation)(factor),
            $message,
            |new_factors| Self::new(new_factors),
            || Ok($crate::expressions::ZeroOperator::new()),
        )
    }};

    // Termwise operation for `Mul`
    (@mul_termwise_operation_all_zero $self:ident, $operation:expr, $message:expr, true) => {{
        $crate::internal::transform_mul_all_zero(
            $self,
            &$self.coefficient.clone().into(),
            &$self.factors,
            |factor: expr_arc_ref_ty!()| ($operation)(factor),
            $message,
            |new_factors| Self::new(new_factors),
            || Ok($crate::expressions::Number::zero()),
        )
    }};

    // Termwise operation for `MatrixMul`
    (@mul_termwise_operation_all_zero $self:ident, $operation:expr, $message:expr, false) => {{
        $crate::internal::transform_mul_all_zero(
            $self,
            &$self.coefficient,
            &$self.factors,
            |factor: expr_arc_ref_ty!()| ($operation)(factor),
            $message,
            |new_factors| Self::new(new_factors),
            || Ok($crate::expressions::ZeroOperator::new()),
        )
    }};

    (@mul_wrote_coefficient $coefficient:expr, true) => {
        !$coefficient.is_one(None)
    };

    (@mul_wrote_coefficient $coefficient:expr, false) => {
        !$crate::public::is_one_expr(&$coefficient, None)
    };
}
