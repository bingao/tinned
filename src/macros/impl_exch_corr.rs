macro_rules! impl_exch_corr_type {
    (
        $type_name:ident,       // ExchCorrEnergy or ExchCorrPotential
        $builder_name:ident,    // ExchCorrEnergyBuilder or ExchCorrPotentialBuilder
        $grid_expr_name:ident,  // xc_energy or xc_potential
        $is_scalar:tt
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

                impl_exch_corr_type!(@build_grid_expr self $grid_expr_name $is_scalar);

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

    (@build_grid_expr $self:ident $grid_expr_name:ident true) => {
        let xc_density = build_xc_density(
            "Exc",
            $self.density_matrix.clone(),
            $self.overlap_distribution.clone(),
            0,
        )?;

        // - `Mul` for unperturbed or the first-order perturbed cases
        // - `Add` for higher-order perturbed case
        let $grid_expr_name = crate::expressions::Mul::new(vec![
            $self.grid_weight.clone(),
            xc_density],
        )?;
    };

    (@build_grid_expr $self:ident $grid_expr_name:ident false) => {
        let xc_density = build_xc_density(
            "Vxc",
            $self.density_matrix.clone(),
            $self.overlap_distribution.clone(),
            1,
        )?;

        // - `MatrixMul` for unperturbed case, or when the generalized overlap
        //   distribution does not depend on applied perturbation(s)
        // - `MatrixAdd` for perturbed case in particular the generalized
        //   overlap distribution depends on applied perturbation(s)
        let $grid_expr_name = crate::expressions::MatrixMul::new(vec![
            crate::expressions::Mul::new(vec![$self.grid_weight.clone(), xc_density])?,
            $self.overlap_distribution.clone()
        ])?;
    };
}

macro_rules! impl_exch_corr_traits {
    ($type_name:ident, $grid_expr_name:ident, $is_scalar:literal) => {
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
            fn eq_expr(&self, other: &dyn Expr) -> bool {
                if let Some(xc) = downcast_from_ref::<$type_name>(other) {
                    self == xc
                } else {
                    false
                }
            }

            #[inline]
            fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "{self}")
            }

            fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
                let diff_expr = self.$grid_expr_name.differentiate(s)?;
                let mut new_deriv = self.derivative.clone();
                new_deriv.insert(s);

                Ok(intern_expr(Arc::new(Self {
                    name: self.name.clone(),
                    grid_weight: self.grid_weight.clone(),
                    density_matrix: self.density_matrix.clone(),
                    overlap_distribution: self.overlap_distribution.clone(),
                    $grid_expr_name: diff_expr,
                    derivative: new_deriv,
                })))
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
