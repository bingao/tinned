macro_rules! impl_exch_corr_type {
    (
        $type_name:ident,          // ExchCorrEnergy or ExchCorrPotential
        $builder_name:ident,       // ExchCorrEnergyBuilder or ExchCorrPotentialBuilder
        $xc_grid_expr_name:ident,  // xc_energy or xc_potential
        $is_scalar:literal
    ) => {
        #[derive(Clone, Debug, Serialize, Deserialize)]
        pub struct $type_name {
            name: String,
            grid_weight: Arc<dyn Expr>,
            density_matrix: Arc<dyn Expr>,
            overlap_distribution: Arc<dyn Expr>,
            $xc_grid_expr_name: Arc<dyn Expr>,
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
                $builder_name { name: name.into(), grid_weight, density_matrix, overlap_distribution }
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
            pub fn $xc_grid_expr_name(&self) -> &Arc<dyn Expr> {
                &self.$xc_grid_expr_name
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

                let order = if $is_scalar { 0 } else { 1 };
                let xc_density = build_xc_density(
                    if $is_scalar { "Exc" } else { "Vxc" },
                    self.density_matrix.clone(),
                    self.overlap_distribution.clone(),
                    order,
                )?;

                let grid_expr = if $is_scalar {
                    Mul::new(vec![self.grid_weight.clone(), xc_density])?
                } else {
                    crate::expressions::MatrixMul::new(vec![
                        Mul::new(vec![self.grid_weight.clone(), xc_density])?,
                        self.overlap_distribution.clone(),
                    ])?
                };

                Ok(intern(Arc::new($type_name {
                    name: self.name,
                    grid_weight: self.grid_weight,
                    density_matrix: self.density_matrix,
                    overlap_distribution: self.overlap_distribution,
                    $xc_grid_expr_name: grid_expr,
                    derivative: PertMultichain::new(),
                })))
            }
        }
    };
}

macro_rules! impl_exch_corr_traits {
    (
        $type_name:ident,          // ExchCorrEnergy or ExchCorrPotential
        $xc_grid_expr_name:ident,  // xc_energy or xc_potential
        $xc_grid_term_type:ty,     // Mul or MatrixMul
        $xc_grid_expr_type:ty,     // Add or MatrixAdd
        $fmt_xc_grid_term:ident,   // fmt_mul or fmt_matrix_mul
        $is_scalar:literal
    ) => {
        impl Expr for $type_name {
            #[inline]
            fn as_any(&self) -> &dyn Any {
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
                    pert_multichain_hash_key(&self.derivative),
                    self.$xc_grid_expr_name.hash_key(),
                )
            }

            #[inline]
            fn is_scalar(&self) -> bool {
                $is_scalar
            }

            #[inline]
            fn eq_expr(&self, other: &dyn Expr) -> bool {
                if let Some(xc) = downcast_expr::<$type_name>(other) {
                    self == xc
                } else {
                    false
                }
            }

            fn differentiate(&self, s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError> {
                let diff_expr = self.$xc_grid_expr_name.differentiate(s)?;
                let mut new_deriv = self.derivative.clone();
                *new_deriv.entry(s.clone()).or_insert(0) += 1;

                Ok(intern(Arc::new(Self {
                    name: self.name.clone(),
                    grid_weight: self.grid_weight.clone(),
                    density_matrix: self.density_matrix.clone(),
                    overlap_distribution: self.overlap_distribution.clone(),
                    $xc_grid_expr_name: diff_expr,
                    derivative: new_deriv,
                })))
            }
        }

        impl PartialEq for $type_name {
            fn eq(&self, other: &Self) -> bool {
                if self.name != other.name
                    || self.grid_weight != other.grid_weight
                    || self.density_matrix != other.density_matrix
                    || self.overlap_distribution != other.overlap_distribution
                    || self.derivative != other.derivative
                {
                    return false;
                }

                self.$xc_grid_expr_name == other.$xc_grid_expr_name
            }
        }

        impl Eq for $type_name {}

        impl Display for $type_name {
            fn fmt(&self, f: &mut Formatter) -> FmtResult {
                // ExchCorrEnergy: unperturbed or the first-order perturbed cases
                //
                // ExchCorrPotential: unperturbed case or when the generalized
                // overlap distribution does not depend on applied
                // perturbation(s)
                if let Some(mul) = downcast_expr::<$xc_grid_term_type>(&self.$xc_grid_expr_name) {
                    write!(f, "{}[", self.name)?;
                    $fmt_xc_grid_term(f, mul)?;
                    write!(f, "]")
                // ExchCorrEnergy: Higher-order perturbed case
                //
                // ExchCorrPotential: perturbed case in particular the
                // generalized overlap distribution depends on applied
                // perturbation(s)
                } else if let Some(add) = downcast_expr::<$xc_grid_expr_type>(&self.$xc_grid_expr_name)
                {
                    let mut first_term = true;
                    for term in add.terms() {
                        if let Some(mul) = downcast_expr::<$xc_grid_term_type>(term) {
                            if !first_term {
                                write!(f, " + ")?;
                            }
                            write!(f, "{}[", self.name)?;
                            $fmt_xc_grid_term(f, mul)?;
                            write!(f, "]")?;
                            first_term = false;
                        } else {
                            write!(f, "{}[unreachable term type: {}]", self.name, term)?;
                        }
                    }
                    Ok(())
                } else {
                    write!(f, "{}[unreachable expr type: {}]", self.name, self.$xc_grid_expr_name)
                }
            }
        }
    };
}
