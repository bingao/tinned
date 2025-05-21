macro_rules! impl_exch_corr_type {
    (
        $type_name:ident,       // ExchCorrEnergy or ExchCorrPotential
        $builder_name:ident,    // ExchCorrEnergyBuilder or ExchCorrPotentialBuilder
        $grid_expr_name:ident,  // xc_energy or xc_potential
        $build_grid_expr:ident  // build_xc_energy or build_xc_potential
    ) => {
        #[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
        pub struct $type_name {
            name: String,
            grid_weight: Arc<dyn Expr>,
            density_matrix: Arc<dyn Expr>,
            overlap_distribution: Arc<dyn Expr>,
            $grid_expr_name: Arc<dyn Expr>,
            derivative: PertMultichain,
        }

        impl $type_name {
            #[inline]
            pub fn builder(
                name: impl Into<String>,
                grid_weight: Arc<dyn Expr>,
                density_matrix: Arc<dyn Expr>,
                overlap_distribution: Arc<dyn Expr>,
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
            pub fn grid_weight(&self) -> &Arc<dyn Expr> {
                &self.grid_weight
            }

            #[inline]
            pub fn density_matrix(&self) -> &Arc<dyn Expr> {
                &self.density_matrix
            }

            #[inline]
            pub fn overlap_distribution(&self) -> &Arc<dyn Expr> {
                &self.overlap_distribution
            }

            #[inline]
            pub fn $grid_expr_name(&self) -> &Arc<dyn Expr> {
                &self.$grid_expr_name
            }

            #[inline]
            pub fn derivative(&self) -> &PertMultichain {
                &self.derivative
            }
        }

        #[derive(Debug)]
        pub struct $builder_name {
            name: String,
            grid_weight: Arc<dyn Expr>,
            density_matrix: Arc<dyn Expr>,
            overlap_distribution: Arc<dyn Expr>,
        }

        impl $builder_name {
            pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
                validate_xc_inputs(
                    &self.density_matrix,
                    &self.grid_weight,
                    &self.overlap_distribution,
                )?;

                let $grid_expr_name = $build_grid_expr(
                    self.grid_weight.clone(),
                    self.density_matrix.clone(),
                    self.overlap_distribution.clone(),
                )?;

                Ok(intern_expr(Arc::new($type_name {
                    name: self.name,
                    grid_weight: self.grid_weight,
                    density_matrix: self.density_matrix,
                    overlap_distribution: self.overlap_distribution,
                    $grid_expr_name,
                    derivative: PertMultichain::new(),
                })))
            }
        }
    };
}

macro_rules! impl_exch_corr_traits {
    ($type_name:ident, $grid_expr_name:ident, $is_scalar:tt) => {
        impl ExprInternal for  $type_name {
            impl_expr_internal_methods!($type_name);

            #[inline]
            fn hash_key(&self) -> String {
                format!(
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
            fn match_for_find_all(&self, other: &Arc<dyn Expr>) -> bool {
                if let Some(xc) = downcast_from_arc::<$type_name>(other) {
                    self.name == xc.name
                        && &self.grid_weight == &xc.grid_weight
                        && &self.density_matrix == &xc.density_matrix
                        && &self.overlap_distribution == &xc.overlap_distribution
                } else {
                    false
                }
            }

            #[inline]
            fn match_for_replace_all(&self, other: &Arc<dyn Expr>) -> bool {
                self.match_for_find_all(other)
            }
        }

        #[typetag::serde]
        impl Expr for $type_name {
            impl_expr_common_methods!($is_scalar);

            fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
                let diff_expr = self.$grid_expr_name.differentiate(s).map_err(|e| {
                    generic_expression_error(
                        concat!(
                            stringify!($type_name),
                            "::differentiate() failed for ",
                            stringify!($grid_expr_name)
                        ),
                        self,
                        Some(Box::new(e)),
                    )
                })?;

                let new_deriv = self.derivative.with_added_perturbation(s);

                Ok(intern_expr(Arc::new(Self {
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
                parameter: &Arc<dyn Expr>,
                perturbations: &[Arc<Perturbation>],
                min_order: u32,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                impl_exch_corr_traits!(
                    @grid_expr_operation
                    self,
                    $grid_expr_name,
                    |grid_expr: &Arc<dyn Expr>| grid_expr.eliminate(parameter, perturbations, min_order),
                    concat!(stringify!($type_name), "::eliminate() failed"),
                    $is_scalar,
                    false
                )
            }

            #[inline]
            fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
                set.iter().any(|expr| self.eq_expr(expr.as_ref()))
                    || self.$grid_expr_name.exist_any(set)
            }

            #[inline]
            fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
                if self.match_for_find_all(s) {
                    BTreeMap::from([(self.total_order(), HashSet::from([self.clone_expr()]))])
                } else {
                    self.$grid_expr_name.find_all(s)
                }
            }

            #[inline]
            fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
                if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                    return impl_zero_expr!($is_scalar);
                }

                impl_exch_corr_traits!(
                    @grid_expr_operation
                    self,
                    $grid_expr_name,
                    |grid_expr: &Arc<dyn Expr>| grid_expr.remove(set),
                    concat!(stringify!($type_name), "::remove() failed"),
                    $is_scalar,
                    false
                )
            }

            #[inline]
            fn retain(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
                if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                    return Ok(self.clone_expr());
                }

                impl_exch_corr_traits!(
                    @grid_expr_operation
                    self,
                    $grid_expr_name,
                    |grid_expr: &Arc<dyn Expr>| grid_expr.retain(set),
                    concat!(stringify!($type_name), "::retain() failed"),
                    $is_scalar,
                    false
                )
            }

            #[inline]
            fn replace(
                &self,
                map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                if let Some((_, value)) = map.iter().find(|(key, _)| self.eq_expr(key.as_ref())) {
                    return if self.derivative.is_empty() {
                        Ok(value.clone())
                    } else {
                        differentiate_expr(value, &self.derivative)
                    };
                }

                impl_exch_corr_traits!(
                    @grid_expr_operation
                    self,
                    $grid_expr_name,
                    |grid_expr: &Arc<dyn Expr>| grid_expr.replace(map),
                    concat!(stringify!($type_name), "::replace() failed"),
                    $is_scalar,
                    true
                )
            }

            #[inline]
            fn replace_all(
                &self,
                map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
            ) -> Result<Arc<dyn Expr>, TinnedError> {
                if let Some((_, value)) = map.iter().find(|(key, _)| self.match_for_replace_all(key)) {
                    return if self.derivative.is_empty() {
                        Ok(value.clone())
                    } else {
                        differentiate_expr(value, &self.derivative)
                    };
                }

                impl_exch_corr_traits!(
                    @grid_expr_operation
                    self,
                    $grid_expr_name,
                    |grid_expr: &Arc<dyn Expr>| grid_expr.replace_all(map),
                    concat!(stringify!($type_name), "::replace_all() failed"),
                    $is_scalar,
                    true
                )
            }
        }

        impl PartialEq for $type_name {
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

        impl Eq for $type_name {}

        impl std::fmt::Display for $type_name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "{}[{}]", self.name, &self.$grid_expr_name)
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
            generic_expression_error(
                concat!($message, " for grid weight"),
                $self,
                Some(Box::new(e))
            )
        })?;
        let density_matrix = ($operation)(&$self.density_matrix).map_err(|e| {
            generic_expression_error(
                concat!($message, " for density matrix"),
                $self,
                Some(Box::new(e))
            )
        })?;
        let overlap_distribution = ($operation)(&$self.overlap_distribution).map_err(|e| {
            generic_expression_error(
                concat!($message, " for overlap distribution"),
                $self,
                Some(Box::new(e))
            )
        })?;

        let new_expr = ($operation)(&$self.$grid_expr_name).map_err(|e| {
            generic_expression_error(
                concat!($message, " for ", stringify!($grid_expr_name)),
                $self,
                Some(Box::new(e)),
            )
        })?;

        if is_zero_expr(&new_expr, None) {
            impl_zero_expr!($is_scalar)
        } else if &new_expr == &$self.$grid_expr_name {
            Ok($self.clone_expr())
        } else {
            Ok(intern_expr(Arc::new(Self {
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
            generic_expression_error(
                concat!($message, " for ", stringify!($grid_expr_name)),
                $self,
                Some(Box::new(e)),
            )
        })?;

        if is_zero_expr(&new_expr, None) {
            impl_zero_expr!($is_scalar)
        } else if &new_expr == &$self.$grid_expr_name {
            Ok($self.clone_expr())
        } else {
            Ok(intern_expr(Arc::new(Self {
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
            grid_weight: Option<Arc<dyn Expr>>,
            density_matrix: Option<Arc<dyn Expr>>,
            overlap_distribution: Option<Arc<dyn Expr>>,
        ) -> Arc<dyn Expr> {
            let name: String = name.into();
            let weight = grid_weight.unwrap_or_else(|| make_non_elec_function(""));
            let dens = density_matrix.unwrap_or_else(|| make_wfn_parameter(""));
            let overlap = overlap_distribution.unwrap_or_else(|| make_one_elec_operator(""));
            if name.is_empty() {
                $type_name::builder(
                    random_alphanumeric($oper_name.len() as u32 + 1),
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
