macro_rules! impl_exch_corr_type {
    (
        $type_name:ident,       // ExchCorrEnergy or ExchCorrPotential
        $builder_name:ident,    // ExchCorrEnergyBuilder or ExchCorrPotentialBuilder
        $grid_expr_name:ident,  // xc_energy or xc_potential
        $build_grid_expr:ident  // build_xc_energy or build_xc_potential
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
        }

        #[derive(Debug)]
        pub struct $builder_name {
            name: String,
            grid_weight: expr_arc_ty!(),
            density_matrix: expr_arc_ty!(),
            overlap_distribution: expr_arc_ty!(),
        }

        impl $builder_name {
            pub fn build(self) -> expr_result_ty!() {
                $crate::internal::validate_xc_inputs(
                    &self.density_matrix,
                    &self.grid_weight,
                    &self.overlap_distribution,
                )?;

                let $grid_expr_name = $build_grid_expr(
                    self.grid_weight.clone(),
                    self.density_matrix.clone(),
                    self.overlap_distribution.clone(),
                )?;

                Ok($crate::internal::intern_expr(::std::sync::Arc::new($type_name {
                    name: self.name,
                    grid_weight: self.grid_weight,
                    density_matrix: self.density_matrix,
                    overlap_distribution: self.overlap_distribution,
                    $grid_expr_name,
                    derivative: $crate::perturbations::PertMultichain::new(),
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
            fn total_order(&self) -> u32 {
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
            fn replace_expr_fields(
                &self,
                map: &expr_map_ty!(),
                exact_equality: bool,
            ) -> expr_result_ty!() {
                impl_exch_corr_traits!(
                    @grid_expr_operation
                    self,
                    $grid_expr_name,
                    |grid_expr: expr_arc_ref_ty!()| grid_expr.replace(map, exact_equality),
                    concat!(stringify!($type_name), "::replace_expr_fields() failed"),
                    $is_scalar,
                    true
                )
            }

            #[inline]
            fn retain_expr_fields(
                &self,
                expr: expr_arc_ref_ty!(),
                exact_equality: bool,
            ) -> expr_result_ty!() {
                impl_exch_corr_traits!(
                    @grid_expr_operation
                    self,
                    $grid_expr_name,
                    |grid_expr: expr_arc_ref_ty!()| grid_expr.retain_expr(expr, exact_equality),
                    concat!(stringify!($type_name), "::retain_expr_fields() failed"),
                    $is_scalar,
                    false
                )
            }
        }

        #[::typetag::serde]
        impl $crate::core::Expr for $type_name {
            impl_expr_common_methods!($is_scalar);

            fn differentiate(
                &self,
                s: &pert_arc_ty!(),
            ) -> expr_result_ty!() {
                let diff_expr = self.$grid_expr_name.differentiate(s).map_err(|e| {
                    $crate::public::generic_expression_error(
                        concat!(
                            stringify!($type_name),
                            "::differentiate() failed for ",
                            stringify!($grid_expr_name)
                        ),
                        self,
                        Some(::std::boxed::Box::new(e)),
                    )
                })?;

                let new_deriv = self.derivative.with_added_perturbation(s);

                Ok($crate::internal::intern_expr(::std::sync::Arc::new(Self {
                    name: self.name.clone(),
                    grid_weight: self.grid_weight.clone(),
                    density_matrix: self.density_matrix.clone(),
                    overlap_distribution: self.overlap_distribution.clone(),
                    $grid_expr_name: diff_expr,
                    derivative: new_deriv,
                })))
            }

            fn eliminate(
                &self,
                parameter: expr_arc_ref_ty!(),
                perturbations: &[pert_arc_ty!()],
                min_order: u32,
            ) -> expr_result_ty!() {
                impl_exch_corr_traits!(
                    @grid_expr_operation
                    self,
                    $grid_expr_name,
                    |grid_expr: expr_arc_ref_ty!()| {
                        grid_expr.eliminate(parameter, perturbations, min_order)
                    },
                    concat!(stringify!($type_name), "::eliminate() failed"),
                    $is_scalar,
                    false
                )
            }

            #[inline]
            fn exist_any(&self, set: &expr_set_ty!()) -> bool {
                set.iter().any(|expr| self.eq_expr(expr.as_ref()))
                    || self.$grid_expr_name.exist_any(set)
            }

            #[inline]
            fn find_superchains(&self, s: expr_arc_ref_ty!()) -> expr_differentiation_map_ty!() {
                if self.deep_eq_superchains(s) {
                    ::std::collections::BTreeMap::from([(
                        self.total_order(),
                        ::std::collections::HashSet::from([self.clone_expr()]),
                    )])
                } else {
                    self.$grid_expr_name.find_superchains(s)
                }
            }

            #[inline]
            fn remove(&self, set: &expr_set_ty!()) -> expr_result_ty!() {
                if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                    return impl_zero_expr!($is_scalar);
                }

                impl_exch_corr_traits!(
                    @grid_expr_operation
                    self,
                    $grid_expr_name,
                    |grid_expr: expr_arc_ref_ty!()| grid_expr.remove(set),
                    concat!(stringify!($type_name), "::remove() failed"),
                    $is_scalar,
                    false
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

    (
        @grid_expr_operation
        $self:ident,
        $grid_expr_name:ident,
        $operation:expr,
        $message:expr,
        $is_scalar:tt,
        true
    ) => {{
        let grid_weight = ($operation)(&$self.grid_weight).map_err(|e| {
            $crate::public::generic_expression_error(
                concat!($message, " for grid weight"),
                $self,
                Some(::std::boxed::Box::new(e)),
            )
        })?;
        let density_matrix = ($operation)(&$self.density_matrix).map_err(|e| {
            $crate::public::generic_expression_error(
                concat!($message, " for density matrix"),
                $self,
                Some(::std::boxed::Box::new(e)),
            )
        })?;
        let overlap_distribution = ($operation)(&$self.overlap_distribution).map_err(|e| {
            $crate::public::generic_expression_error(
                concat!($message, " for overlap distribution"),
                $self,
                Some(::std::boxed::Box::new(e)),
            )
        })?;

        let new_expr = ($operation)(&$self.$grid_expr_name).map_err(|e| {
            $crate::public::generic_expression_error(
                concat!($message, " for ", stringify!($grid_expr_name)),
                $self,
                Some(::std::boxed::Box::new(e)),
            )
        })?;

        if $crate::public::is_zero_expr(&new_expr, None) {
            impl_zero_expr!($is_scalar)
        } else if &new_expr == &$self.$grid_expr_name {
            Ok($self.clone_expr())
        } else {
            Ok($crate::internal::intern_expr(::std::sync::Arc::new(Self {
                name: $self.name.clone(),
                grid_weight,
                density_matrix,
                overlap_distribution,
                $grid_expr_name: new_expr,
                derivative: $self.derivative.clone(),
            })))
        }
    }};

    (
        @grid_expr_operation
        $self:ident,
        $grid_expr_name:ident,
        $operation:expr,
        $message:expr,
        $is_scalar:tt,
        false
    ) => {{
        let new_expr = ($operation)(&$self.$grid_expr_name).map_err(|e| {
            $crate::public::generic_expression_error(
                concat!($message, " for ", stringify!($grid_expr_name)),
                $self,
                Some(::std::boxed::Box::new(e)),
            )
        })?;

        if $crate::public::is_zero_expr(&new_expr, None) {
            impl_zero_expr!($is_scalar)
        } else if &new_expr == &$self.$grid_expr_name {
            Ok($self.clone_expr())
        } else {
            Ok($crate::internal::intern_expr(::std::sync::Arc::new(Self {
                name: $self.name.clone(),
                grid_weight: $self.grid_weight.clone(),
                density_matrix: $self.density_matrix.clone(),
                overlap_distribution: $self.overlap_distribution.clone(),
                $grid_expr_name: new_expr,
                derivative: $self.derivative.clone(),
            })))
        }
    }};
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
                $crate::expressions::one_elec_operator::test_utils::make_one_elec_operator("")
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
