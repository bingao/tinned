macro_rules! impl_exch_corr_type {
    (
        $type_name:ident,       // ExchCorrEnergy or ExchCorrPotential
        $builder_name:ident,    // ExchCorrEnergyBuilder or ExchCorrPotentialBuilder
        $grid_expr_name:ident,  // xc_energy or xc_potential
        $build_grid_expr:ident, // build_xc_energy or build_xc_potential
        $is_scalar:tt           // Whether the expression is scalar or not
    ) => {
        #[derive(Clone, Debug, ::serde::Serialize, ::serde::Deserialize)]
        pub struct $type_name {
            name: String,
            grid_weight: expr_arc_ty!(),
            density_matrix: expr_arc_ty!(),
            overlap_distribution: expr_arc_ty!(),
            $grid_expr_name: expr_arc_ty!(),
            derivative: $crate::perturbations::PertMultichain,
        }

        impl $type_name {
            #[inline]
            pub fn builder(
                name: impl Into<String>,
                grid_weight: expr_arc_ty!(),
                density_matrix: expr_arc_ty!(),
                overlap_distribution: expr_arc_ty!(),
            ) -> $builder_name {
                $builder_name {
                    name: name.into(),
                    grid_weight,
                    density_matrix,
                    overlap_distribution,
                    $grid_expr_name: None,
                    derivative: None,
                }
            }

            #[inline]
            fn with_grid_expr(&self, grid_expr: expr_arc_ty!()) -> $builder_name {
                $builder_name {
                    name: self.name.clone(),
                    grid_weight: self.grid_weight.clone(),
                    density_matrix: self.density_matrix.clone(),
                    overlap_distribution: self.overlap_distribution.clone(),
                    $grid_expr_name: Some(grid_expr),
                    derivative: Some(self.derivative.clone()),
                }
            }

            #[inline]
            fn with_grid_expr_and_derivative(
                &self,
                grid_expr: expr_arc_ty!(),
                derivative: $crate::perturbations::PertMultichain,
            ) -> $builder_name {
                $builder_name {
                    name: self.name.clone(),
                    grid_weight: self.grid_weight.clone(),
                    density_matrix: self.density_matrix.clone(),
                    overlap_distribution: self.overlap_distribution.clone(),
                    $grid_expr_name: Some(grid_expr),
                    derivative: Some(derivative),
                }
            }

            #[inline]
            pub fn name(&self) -> &str {
                &self.name
            }

            #[inline]
            pub fn grid_weight(&self) -> expr_arc_ref_ty!() {
                &self.grid_weight
            }

            #[inline]
            pub fn density_matrix(&self) -> expr_arc_ref_ty!() {
                &self.density_matrix
            }

            #[inline]
            pub fn overlap_distribution(&self) -> expr_arc_ref_ty!() {
                &self.overlap_distribution
            }

            #[inline]
            pub fn $grid_expr_name(&self) -> expr_arc_ref_ty!() {
                &self.$grid_expr_name
            }

            #[inline]
            pub fn derivative(&self) -> &$crate::perturbations::PertMultichain {
                &self.derivative
            }

            fn apply_to_children(
                &self,
                operation: impl Fn(&Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
                message: impl Into<String>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                let message = message.into();

                let grid_weight = operation(&self.grid_weight).map_err(|e| {
                    $crate::public::generic_expression_error(
                        format!("{} for grid weight", message),
                        self,
                        Some(::std::boxed::Box::new(e)),
                    )
                })?;

                let density_matrix = operation(&self.density_matrix).map_err(|e| {
                    $crate::public::generic_expression_error(
                        format!("{} for density matrix", message),
                        self,
                        Some(::std::boxed::Box::new(e)),
                    )
                })?;

                let overlap_distribution = operation(&self.overlap_distribution).map_err(|e| {
                    $crate::public::generic_expression_error(
                        format!("{} for overlap distribution", message),
                        self,
                        Some(::std::boxed::Box::new(e)),
                    )
                })?;

                $crate::internal::transform_unary_any_zero(
                    self,
                    &self.$grid_expr_name,
                    |grid_expr| operation(grid_expr),
                    format!("{} for {}", message, stringify!($grid_expr_name)),
                    |grid_expr| {
                        $type_name::builder(
                            self.name.clone(),
                            grid_weight,
                            density_matrix,
                            overlap_distribution,
                        )
                        .$grid_expr_name(grid_expr)
                        .derivative(self.derivative.clone())
                        .build()
                    },
                    || impl_zero_expr!($is_scalar),
                )
            }
        }

        #[derive(Debug)]
        pub struct $builder_name {
            name: String,
            grid_weight: expr_arc_ty!(),
            density_matrix: expr_arc_ty!(),
            overlap_distribution: expr_arc_ty!(),
            $grid_expr_name: ::std::option::Option<expr_arc_ty!()>,
            derivative: ::std::option::Option<$crate::perturbations::PertMultichain>,
        }

        impl $builder_name {
            #[inline]
            fn $grid_expr_name(mut self, $grid_expr_name: expr_arc_ty!()) -> Self {
                self.$grid_expr_name = Some($grid_expr_name);
                self
            }

            #[inline]
            fn derivative(mut self, derivative: $crate::perturbations::PertMultichain) -> Self {
                self.derivative = Some(derivative);
                self
            }

            pub fn build(self) -> expr_result_ty!() {
                $crate::internal::validate_xc_inputs(
                    &self.density_matrix,
                    &self.grid_weight,
                    &self.overlap_distribution,
                )?;

                let $grid_expr_name = self.$grid_expr_name.unwrap_or($build_grid_expr(
                    self.grid_weight.clone(),
                    self.density_matrix.clone(),
                    self.overlap_distribution.clone(),
                )?);

                let derivative =
                    self.derivative.unwrap_or($crate::perturbations::PertMultichain::new());

                Ok($crate::internal::intern_expr(::std::sync::Arc::new($type_name {
                    name: self.name,
                    grid_weight: self.grid_weight,
                    density_matrix: self.density_matrix,
                    overlap_distribution: self.overlap_distribution,
                    $grid_expr_name,
                    derivative,
                })))
            }
        }
    };
}

macro_rules! impl_exch_corr_traits {
    ($type_name:ident, $grid_expr_name:ident, $is_scalar:tt) => {
        impl $crate::core::ExprInternal for $type_name {
            impl_expr_internal_methods!($type_name, true);

            #[inline]
            fn hash_key(&self) -> String {
                ::std::format!(
                    "{}({}; {}; {}; {}; [{}]; {})",
                    stringify!($type_name),
                    self.name,
                    self.grid_weight.hash_key(),
                    self.density_matrix.hash_key(),
                    self.overlap_distribution.hash_key(),
                    self.derivative.hash_key(),
                    self.$grid_expr_name.hash_key(),
                )
            }

            #[inline]
            fn expr_order(&self) -> u32 {
                self.derivative.total_order()
            }

            #[inline]
            fn deep_eq_superchains(&self, other: expr_arc_ref_ty!()) -> bool {
                if let Some(xc) = $crate::public::downcast_from_arc::<$type_name>(other) {
                    self.name == xc.name
                        && &self.grid_weight == &xc.grid_weight
                        && &self.density_matrix == &xc.density_matrix
                        && &self.overlap_distribution == &xc.overlap_distribution
                        && self.derivative.is_subchain(&xc.derivative)
                } else {
                    false
                }
            }

            #[inline]
            fn eq_by_superchains(&self, other: expr_arc_ref_ty!()) -> bool {
                self.deep_eq_superchains(other)
            }

            #[inline]
            fn replace_one_in_children(
                &self,
                expr: expr_arc_ref_ty!(),
                replacement: expr_arc_ty!(),
                include_derivatives: bool,
            ) -> expr_result_ty!() {
                self.apply_to_children(
                    |arg: expr_arc_ref_ty!()| {
                        arg.replace_one(expr, replacement.clone(), include_derivatives)
                    },
                    concat!(stringify!($type_name), "::replace_one_in_children() failed"),
                )
            }

            #[inline]
            fn replace_all_in_children(
                &self,
                map: &expr_map_ty!(),
                include_derivatives: bool,
            ) -> expr_result_ty!() {
                self.apply_to_children(
                    |arg: expr_arc_ref_ty!()| arg.replace_all(map, include_derivatives),
                    concat!(stringify!($type_name), "::replace_all_in_children() failed"),
                )
            }
        }

        #[::typetag::serde]
        impl $crate::core::Expr for $type_name {
            impl_expr_common_methods!($is_scalar);

            #[inline]
            fn has_unperturbed_term(&self) -> bool {
                self.grid_weight.has_unperturbed_term()
                    && self.density_matrix.has_unperturbed_term()
                    && self.overlap_distribution.has_unperturbed_term()
            }

            #[inline]
            fn substitute_zero_perturbations(
                &self,
                freq_tol: ::std::option::Option<$crate::public::NumberTolerance>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                $crate::internal::transform_unary_any_zero(
                    self,
                    &self.$grid_expr_name,
                    |grid_expr| grid_expr.substitute_zero_perturbations(freq_tol),
                    concat!(
                        stringify!($type_name),
                        "::substitute_zero_perturbations() failed for ",
                        stringify!($grid_expr_name)
                    ),
                    |grid_expr| self.with_grid_expr(grid_expr).build(),
                    || impl_zero_expr!($is_scalar),
                )
            }

            fn differentiate(&self, s: pert_arc_ty!()) -> expr_result_ty!() {
                let derivative = self.derivative.with_added_perturbation(s.clone());

                $crate::internal::transform_unary_any_zero(
                    self,
                    &self.$grid_expr_name,
                    |grid_expr| grid_expr.differentiate(s),
                    concat!(
                        stringify!($type_name),
                        "::differentiate() failed for ",
                        stringify!($grid_expr_name)
                    ),
                    |grid_expr| self.with_grid_expr_and_derivative(grid_expr, derivative).build(),
                    || impl_zero_expr!($is_scalar),
                )
            }

            fn eliminate(
                &self,
                parameter: expr_arc_ty!(),
                perturbations: &[pert_arc_ty!()],
                min_order: u32,
            ) -> expr_result_ty!() {
                $crate::internal::transform_unary_any_zero(
                    self,
                    &self.$grid_expr_name,
                    |grid_expr| grid_expr.eliminate(parameter, perturbations, min_order),
                    concat!(
                        stringify!($type_name),
                        "::eliminate() failed for ",
                        stringify!($grid_expr_name)
                    ),
                    |grid_expr| self.with_grid_expr(grid_expr).build(),
                    || impl_zero_expr!($is_scalar),
                )
            }

            #[inline]
            fn find_all(&self, s: expr_arc_ref_ty!()) -> expr_differentiation_map_ty!() {
                if self.deep_eq_superchains(s) {
                    ::std::collections::BTreeMap::from([(
                        self.expr_order(),
                        ::std::collections::HashSet::from([self.clone_expr()]),
                    )])
                } else {
                    self.$grid_expr_name.find_all(s)
                }
            }

            #[inline]
            fn match_one(&self, s: expr_arc_ref_ty!(), include_derivatives: bool) -> bool {
                self.match_one_self(s, include_derivatives)
                    || self.$grid_expr_name.match_one(s, include_derivatives)
            }

            #[inline]
            fn match_any(&self, set: &expr_set_ty!(), include_derivatives: bool) -> bool {
                self.match_any_self(set, include_derivatives)
                    || self.$grid_expr_name.match_any(set, include_derivatives)
            }

            #[inline]
            fn remove_one(&self, s: expr_arc_ref_ty!()) -> expr_result_ty!() {
                if self.match_one_self(s, false) {
                    return impl_zero_expr!($is_scalar);
                }

                $crate::internal::transform_unary_any_zero(
                    self,
                    &self.$grid_expr_name,
                    |grid_expr| grid_expr.remove_one(s),
                    concat!(
                        stringify!($type_name),
                        "::remove_one() failed for ",
                        stringify!($grid_expr_name)
                    ),
                    |grid_expr| self.with_grid_expr(grid_expr).build(),
                    || impl_zero_expr!($is_scalar),
                )
            }

            #[inline]
            fn remove_all(&self, set: &expr_set_ty!()) -> expr_result_ty!() {
                if self.match_any_self(set, false) {
                    return impl_zero_expr!($is_scalar);
                }

                $crate::internal::transform_unary_any_zero(
                    self,
                    &self.$grid_expr_name,
                    |grid_expr| grid_expr.remove_all(set),
                    concat!(
                        stringify!($type_name),
                        "::remove_all() failed for ",
                        stringify!($grid_expr_name)
                    ),
                    |grid_expr| self.with_grid_expr(grid_expr).build(),
                    || impl_zero_expr!($is_scalar),
                )
            }

            #[inline]
            fn retain_one(
                &self,
                s: expr_arc_ref_ty!(),
                include_derivatives: bool,
            ) -> expr_result_ty!() {
                if self.match_one_self(s, include_derivatives) {
                    return Ok(self.clone_expr());
                }

                $crate::internal::transform_unary_any_zero(
                    self,
                    &self.$grid_expr_name,
                    |grid_expr| grid_expr.retain_one(s, include_derivatives),
                    concat!(
                        stringify!($type_name),
                        "::retain_one() failed for ",
                        stringify!($grid_expr_name)
                    ),
                    |grid_expr| self.with_grid_expr(grid_expr).build(),
                    || impl_zero_expr!($is_scalar),
                )
            }

            #[inline]
            fn retain_any(
                &self,
                set: &expr_set_ty!(),
                include_derivatives: bool,
            ) -> expr_result_ty!() {
                if self.match_any_self(set, include_derivatives) {
                    return Ok(self.clone_expr());
                }

                $crate::internal::transform_unary_any_zero(
                    self,
                    &self.$grid_expr_name,
                    |grid_expr| grid_expr.retain_any(set, include_derivatives),
                    concat!(
                        stringify!($type_name),
                        "::retain_any() failed for ",
                        stringify!($grid_expr_name)
                    ),
                    |grid_expr| self.with_grid_expr(grid_expr).build(),
                    || impl_zero_expr!($is_scalar),
                )
            }
        }

        impl ::std::cmp::PartialEq for $type_name {
            fn eq(&self, other: &Self) -> bool {
                if self.name != other.name
                    || &self.grid_weight != &other.grid_weight
                    || &self.density_matrix != &other.density_matrix
                    || &self.overlap_distribution != &other.overlap_distribution
                    || self.derivative != other.derivative
                {
                    return false;
                }

                &self.$grid_expr_name == &other.$grid_expr_name
            }
        }

        impl ::std::cmp::Eq for $type_name {}

        impl ::std::fmt::Display for $type_name {
            fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
                ::std::write!(f, "{}[{}]", self.name, &self.$grid_expr_name)
            }
        }
    };
}

#[allow(unused_macros)]
macro_rules! impl_exch_corr_test_utils {
    (
        $type_name:ident,  // ExchCorrEnergy or ExchCorrPotential
        $oper_name:ident,  // DEFAULT_FUNC_NAME
        $make_expr:ident   // make_exch_corr_energy or make_exch_corr_potential
    ) => {
        #[inline]
        pub fn $make_expr(
            name: impl Into<String>,
            grid_weight: ::std::option::Option<expr_arc_ty!()>,
            density_matrix: ::std::option::Option<expr_arc_ty!()>,
            overlap_distribution: ::std::option::Option<expr_arc_ty!()>,
        ) -> expr_arc_ty!() {
            let name: String = name.into();
            let weight = grid_weight.unwrap_or_else(|| {
                $crate::expressions::non_elec_function::test_utils::make_non_elec_function("")
            });
            let dens = density_matrix.unwrap_or_else(|| {
                $crate::expressions::wfn_parameter::test_utils::make_wfn_parameter("")
            });
            let overlap = overlap_distribution.unwrap_or_else(|| {
                $crate::expressions::one_elec_matrix::test_utils::make_one_elec_matrix("", false)
            });

            if name.is_empty() {
                $type_name::builder(
                    $crate::expressions::symbol::test_utils::random_alphanumeric(
                        $oper_name.len() as u32 + 1,
                    ),
                    weight,
                    dens,
                    overlap,
                )
                .build()
                .unwrap()
            } else {
                $type_name::builder(name, weight, dens, overlap).build().unwrap()
            }
        }
    };
}
