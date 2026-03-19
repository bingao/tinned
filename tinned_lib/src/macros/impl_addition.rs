macro_rules! impl_add_traits {
    ($type_name:ident, $hash_delimiter:ident, $fmt_delimiter:ident, $is_scalar:tt) => {
        impl $crate::core::ExprInternal for $type_name {
            impl_expr_internal_methods!($type_name, false);

            #[inline]
            fn hash_key(&self) -> ::std::string::String {
                format!(
                    "{}({})",
                    stringify!($type_name),
                    $crate::internal::join_mapped(
                        &self.terms,
                        $hash_delimiter,
                        |term| term.hash_key(),
                    ),
                )
            }

            // For unambiguous replacement, we require equality for the
            // whole `Add` so that we do not override methods
            // `eq_by_superchains()` and `replace_expr_self()` of
            // `ExprInternal`.
            fn replace_expr_fields(
                &self,
                map: &expr_map_ty!(),
                exact_equality: bool,
            ) -> expr_result_ty!() {
                impl_add_traits!(
                    @add_termwise_operation
                    self,
                    |term: expr_arc_ref_ty!()| term.replace(map, exact_equality),
                    concat!(stringify!($type_name), "::replace_expr_fields() failed"),
                    $is_scalar
                )
            }

            fn retain_expr_fields(
                &self,
                expr: expr_arc_ref_ty!(),
                exact_equality: bool,
            ) -> expr_result_ty!() {
                impl_add_traits!(
                    @add_termwise_operation
                    self,
                    |term: expr_arc_ref_ty!()| term.retain_expr(expr, exact_equality),
                    concat!(stringify!($type_name), "::retain_expr_fields() failed"),
                    $is_scalar
                )
            }
        }

        #[::typetag::serde]
        impl $crate::core::Expr for $type_name {
            impl_expr_common_methods!($is_scalar);

            fn apply_zero_rules(
                &self,
                freq_tol: ::std::option::Option<$crate::public::NumberTolerance>,
            ) -> expr_result_ty!() {
                impl_add_traits!(
                    @add_termwise_operation
                    self,
                    |term: expr_arc_ref_ty!()| term.apply_zero_rules(freq_tol.clone()),
                    concat!(stringify!($type_name), "::apply_zero_rules() failed"),
                    $is_scalar
                )
            }

            fn differentiate(
                &self,
                s: &::std::sync::Arc<$crate::perturbations::Perturbation>,
            ) -> expr_result_ty!() {
                let mut diff_terms = ::std::vec::Vec::with_capacity(self.terms.len());

                for term in &self.terms {
                    let diff = term.differentiate(s).map_err(|e| {
                        $crate::public::generic_expression_error(
                            concat!(stringify!($type_name), "::differentiate() failed"),
                            self,
                            Some(::std::boxed::Box::new(e)),
                        )
                    })?;
                    if !$crate::public::is_zero_expr(&diff, None) {
                        diff_terms.push(diff);
                    }
                }

                Self::new(diff_terms)
            }

            fn eliminate(
                &self,
                parameter: expr_arc_ref_ty!(),
                perturbations: &[::std::sync::Arc<$crate::perturbations::Perturbation>],
                min_order: u32,
            ) -> expr_result_ty!() {
                impl_add_traits!(
                    @add_termwise_operation
                    self,
                    |term: expr_arc_ref_ty!()| term.eliminate(parameter, perturbations, min_order),
                    concat!(stringify!($type_name), "::eliminate() failed"),
                    $is_scalar
                )
            }

            #[inline]
            fn exist_any(&self, set: &expr_set_ty!()) -> bool {
                if self.terms.iter().any(|term| term.exist_any(set)) {
                    return true;
                }

                set.iter().any(|expr| self.eq_expr(expr.as_ref()))
            }

            fn find_superchains(&self, s: expr_arc_ref_ty!()) -> expr_differentiation_map_ty!() {
                if self.deep_eq_superchains(s) {
                    return ::std::collections::BTreeMap::from([(
                        self.total_order(),
                        ::std::collections::HashSet::from([self.clone_expr()]),
                    )]);
                }

                let mut result: expr_differentiation_map_ty!() = ::std::collections::BTreeMap::new();

                for term in &self.terms {
                    for (order, subset) in term.find_superchains(s) {
                        result.entry(order).or_default().extend(subset);
                    }
                }

                result
            }

            fn remove(&self, set: &expr_set_ty!()) -> expr_result_ty!() {
                if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                    return impl_zero_expr!($is_scalar);
                }

                impl_add_traits!(
                    @add_termwise_operation
                    self,
                    |term: expr_arc_ref_ty!()| term.remove(set),
                    concat!(stringify!($type_name), "::remove() failed"),
                    $is_scalar
                )
            }
        }

        impl ::std::cmp::PartialEq for $type_name {
            fn eq(&self, other: &Self) -> bool {
                self.terms == other.terms
            }
        }

        impl ::std::cmp::Eq for $type_name {}

        impl ::std::fmt::Display for $type_name {
            fn fmt(
                &self,
                f: &mut ::std::fmt::Formatter,
            ) -> ::std::fmt::Result {
                ::std::write!(
                    f,
                    "({})",
                    $crate::internal::join_mapped(
                        &self.terms,
                        $fmt_delimiter,
                        |term| term.to_string(),
                    )
                )
            }
        }
    };

    (@add_termwise_operation $self:ident, $operation:expr, $message:expr, $is_scalar:tt) => {{
        let mut new_terms = ::std::vec::Vec::with_capacity($self.terms.len());
        let mut new_add = false;

        for term in &$self.terms {
            let new_term = ($operation)(term).map_err(|e| {
                $crate::public::generic_expression_error(
                    $message,
                    $self,
                    Some(::std::boxed::Box::new(e)),
                )
            })?;
            if $crate::public::is_zero_expr(&new_term, None) {
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
