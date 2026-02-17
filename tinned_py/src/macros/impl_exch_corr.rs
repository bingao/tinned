macro_rules! impl_exch_corr_interface {
    (
        type_name = $type_name:ty,
        type_label = $type_label:literal,
        build_fn = $build_fn:ident,
        name_fn = $name_fn:ident,
        grid_weight_fn = $grid_weight_fn:ident,
        density_matrix_fn = $density_matrix_fn:ident,
        overlap_distribution_fn = $overlap_distribution_fn:ident,
        grid_expr_fn = $grid_expr_fn:ident,
        derivative_fn = $derivative_fn:ident,
        register_fn = $register_fn:ident,
        // method on the Rust type: xc_energy or xc_potential
        grid_expr_method = $grid_expr_method:ident,
        // error message when downcast fails
        downcast_err_msg_prefix = $downcast_err_msg_prefix:literal
    ) => {
        #[doc = concat!(
            "Build an ", $type_label, " expression.\n\n",
            "Args:\n",
            "  name: Name of the XC functional.\n",
            "  grid_weight: Grid weights (expression).\n",
            "  density_matrix: Density matrix (expression).\n",
            "  overlap_distribution: Overlap distribution (expression).\n\n",
            "Returns:\n",
            "  A PyExpr wrapping the constructed expression (interned)."
        )]
        #[pyo3::prelude::pyfunction]
        pub fn $build_fn(
            name: ::std::string::String,
            grid_weight: $crate::core::expr::PyExpr,
            density_matrix: $crate::core::expr::PyExpr,
            overlap_distribution: $crate::core::expr::PyExpr,
        ) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let out = <$type_name>::builder(
                name,
                grid_weight.inner().clone(),
                density_matrix.inner().clone(),
                overlap_distribution.inner().clone(),
            )
            .build()
            .map_err($crate::core::errors::to_pyerr)?;

            Ok($crate::core::expr::PyExpr::new(out))
        }

        #[doc = concat!("Return the name of an ", $type_label, " expression.")]
        #[pyo3::prelude::pyfunction]
        pub fn $name_fn(expr: $crate::core::expr::PyExpr) -> ::pyo3::PyResult<::std::string::String> {
            let inner = expr.inner().clone();

            let xc_ref = inner.as_any().downcast_ref::<$type_name>().ok_or_else(|| {
                $crate::core::errors::to_pyerr(tinned::TinnedError::ExpressionError {
                    message: concat!($downcast_err_msg_prefix, "name() expected an ", $type_label, " expression"),
                    expression: inner.to_string(),
                    source: None,
                })
            })?;

            Ok(xc_ref.name().to_string())
        }

        #[doc = concat!("Return the grid_weight of an ", $type_label, " expression.")]
        #[pyo3::prelude::pyfunction]
        pub fn $grid_weight_fn(expr: $crate::core::expr::PyExpr) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let inner = expr.inner().clone();

            let xc_ref = inner.as_any().downcast_ref::<$type_name>().ok_or_else(|| {
                $crate::core::errors::to_pyerr(tinned::TinnedError::ExpressionError {
                    message: concat!($downcast_err_msg_prefix, "grid_weight() expected an ", $type_label, " expression"),
                    expression: inner.to_string(),
                    source: None,
                })
            })?;

            Ok($crate::core::expr::PyExpr::new(xc_ref.grid_weight().clone()))
        }

        #[doc = concat!("Return the density_matrix of an ", $type_label, " expression.")]
        #[pyo3::prelude::pyfunction]
        pub fn $density_matrix_fn(expr: $crate::core::expr::PyExpr) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let inner = expr.inner().clone();

            let xc_ref = inner.as_any().downcast_ref::<$type_name>().ok_or_else(|| {
                $crate::core::errors::to_pyerr(tinned::TinnedError::ExpressionError {
                    message: concat!($downcast_err_msg_prefix, "density_matrix() expected an ", $type_label, " expression"),
                    expression: inner.to_string(),
                    source: None,
                })
            })?;

            Ok($crate::core::expr::PyExpr::new(xc_ref.density_matrix().clone()))
        }

        #[doc = concat!("Return the overlap_distribution of an ", $type_label, " expression.")]
        #[pyo3::prelude::pyfunction]
        pub fn $overlap_distribution_fn(
            expr: $crate::core::expr::PyExpr,
        ) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let inner = expr.inner().clone();

            let xc_ref = inner.as_any().downcast_ref::<$type_name>().ok_or_else(|| {
                $crate::core::errors::to_pyerr(tinned::TinnedError::ExpressionError {
                    message: concat!($downcast_err_msg_prefix, "overlap_distribution() expected an ", $type_label, " expression"),
                    expression: inner.to_string(),
                    source: None,
                })
            })?;

            Ok($crate::core::expr::PyExpr::new(xc_ref.overlap_distribution().clone()))
        }

        #[doc = concat!("Return the grid expression (", stringify!($grid_expr_method), ") of an ", $type_label, " expression.")]
        #[pyo3::prelude::pyfunction]
        pub fn $grid_expr_fn(expr: $crate::core::expr::PyExpr) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let inner = expr.inner().clone();

            let xc_ref = inner.as_any().downcast_ref::<$type_name>().ok_or_else(|| {
                $crate::core::errors::to_pyerr(tinned::TinnedError::ExpressionError {
                    message: concat!($downcast_err_msg_prefix, "grid_expr() expected an ", $type_label, " expression"),
                    expression: inner.to_string(),
                    source: None,
                })
            })?;

            Ok($crate::core::expr::PyExpr::new(xc_ref.$grid_expr_method().clone()))
        }

        #[doc = concat!("Return the derivative of an ", $type_label, " expression.")]
        #[pyo3::prelude::pyfunction]
        pub fn $derivative_fn(
            expr: $crate::core::expr::PyExpr,
        ) -> ::pyo3::PyResult<$crate::perturbations::pert_multichain::PyPertMultichain> {
            let inner = expr.inner().clone();

            let xc_ref = inner.as_any().downcast_ref::<$type_name>().ok_or_else(|| {
                $crate::core::errors::to_pyerr(tinned::TinnedError::ExpressionError {
                    message: concat!($downcast_err_msg_prefix, "derivative() expected an ", $type_label, " expression"),
                    expression: inner.to_string(),
                    source: None,
                })
            })?;

            Ok($crate::perturbations::pert_multichain::PyPertMultichain::new(
                xc_ref.derivative().clone(),
            ))
        }

        pub fn $register_fn(m: &::pyo3::Bound<'_, ::pyo3::types::PyModule>) -> ::pyo3::PyResult<()> {
            m.add_function(::pyo3::wrap_pyfunction!($build_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($name_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($grid_weight_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($density_matrix_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($overlap_distribution_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($grid_expr_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($derivative_fn, m)?)?;
            Ok(())
        }
    };
}
