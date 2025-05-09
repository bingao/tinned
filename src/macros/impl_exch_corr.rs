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
    ($type_name:ident, $grid_expr_name:ident, $is_scalar:literal, $build_zero_expr:expr) => {
        #[typetag::serde]
        impl Expr for $type_name {
            #[inline]
            fn as_any(&self) -> &dyn std::any::Any {
                self
            }

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
            fn is_scalar(&self) -> bool {
                $is_scalar
            }

            #[inline]
            fn clone_expr(&self) -> Arc<dyn Expr> {
                Arc::new(self.clone())
            }

            #[inline]
            fn eq_expr(&self, other: &dyn Expr) -> bool {
                if let Some(xc) = downcast_from_ref::<$type_name>(other) {
                    self == xc
                } else {
                    false
                }
            }

            #[inline]
            fn eq_shallow(&self, other: &dyn Expr) -> bool {
                if let Some(xc) = downcast_from_ref::<$type_name>(other) {
                    self.name == xc.name
                        && self.grid_weight.eq_shallow(xc.grid_weight.as_ref())
                        && self.density_matrix.eq_shallow(xc.density_matrix.as_ref())
                        && self.overlap_distribution.eq_shallow(xc.overlap_distribution.as_ref())
                } else {
                    false
                }
            }

            #[inline]
            fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "{self}")
            }

            fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
                let diff_expr = self.$grid_expr_name.differentiate(s).map_err(|e| {
                    generic_expression_error(
                        concat!(stringify!($type_name), "differentiate() failed"),
                        self,
                        Some(Box::new(e)),
                    )
                })?;

                let new_deriv = self.derivative.clone_with_insert(s);

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
                let new_expr = self
                    .$grid_expr_name
                    .eliminate(parameter, perturbations, min_order)
                    .map_err(|e| {
                        generic_expression_error(
                            concat!(stringify!($type_name), "eliminate() failed"),
                            self,
                            Some(Box::new(e)),
                        )
                    })?;

                if is_zero_expr(&new_expr, None) {
                    Ok($build_zero_expr())
                } else if &new_expr == &self.$grid_expr_name {
                    Ok(self.clone_expr())
                } else {
                    Ok(intern_expr(Arc::new(Self {
                        name: self.name.clone(),
                        grid_weight: self.grid_weight.clone(),
                        density_matrix: self.density_matrix.clone(),
                        overlap_distribution: self.overlap_distribution.clone(),
                        $grid_expr_name: new_expr,
                        derivative: self.derivative.clone(),
                    })))
                }
            }

            #[inline]
            fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
                set.iter().any(|expr| self.eq_expr(expr.as_ref()))
                    || self.$grid_expr_name.exist_any(set)
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
