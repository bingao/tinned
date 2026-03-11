macro_rules! impl_exch_corr_interface {
    (
        expr_ty = $expr_ty:ty,
        build_fn = $build_fn:ident,
        name_fn = $name_fn:ident,
        grid_weight_fn = $grid_weight_fn:ident,
        density_matrix_fn = $density_matrix_fn:ident,
        overlap_distribution_fn = $overlap_distribution_fn:ident,
        grid_expr_fn = $grid_expr_fn:ident,
        derivative_fn = $derivative_fn:ident,
        register_fn = $register_fn:ident,
        // method on the Rust type: xc_energy or xc_potential
        grid_expr_method = $grid_expr_method:ident
    ) => {
        #[doc = concat!(
                            "Build an ", stringify!($expr_ty), " expression.\n\n",
                            "Args:\n",
                            "  name: Name of the XC functional.\n",
                            "  grid_weight: Grid weights (expression).\n",
                            "  density_matrix: Density matrix (expression).\n",
                            "  overlap_distribution: Overlap distribution (expression).\n\n",
                            "Returns:\n",
                            "  A PyExpr wrapping the constructed expression (interned)."
                        )]
        #[::pyo3::prelude::pyfunction]
        pub fn $build_fn(
            name: ::std::string::String,
            grid_weight: $crate::core::expr::PyExpr,
            density_matrix: $crate::core::expr::PyExpr,
            overlap_distribution: $crate::core::expr::PyExpr,
        ) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let out = <$expr_ty>::builder(
                name,
                grid_weight.inner().clone(),
                density_matrix.inner().clone(),
                overlap_distribution.inner().clone(),
            )
            .build()
            .map_err($crate::core::errors::to_pyerr)?;

            Ok($crate::core::expr::PyExpr::new(out))
        }

        impl_expr_getter_interface!(
            fn_name = $name_fn,
            fn_doc = impl_expr_getter_doc!("name", $expr_ty),
            expr_ty = $expr_ty,
            out_ty = ::std::string::String,
            body = |op: &$expr_ty| Ok(op.name().to_string())
        );

        impl_expr_getter_interface!(
            fn_name = $grid_weight_fn,
            fn_doc = impl_expr_getter_doc!("grid weight", $expr_ty),
            expr_ty = $expr_ty,
            out_ty = $crate::core::expr::PyExpr,
            body = |op: &$expr_ty| Ok($crate::core::expr::PyExpr::new(op.grid_weight().clone()))
        );

        impl_expr_getter_interface!(
            fn_name = $density_matrix_fn,
            fn_doc = impl_expr_getter_doc!("density matrix", $expr_ty),
            expr_ty = $expr_ty,
            out_ty = $crate::core::expr::PyExpr,
            body = |op: &$expr_ty| Ok($crate::core::expr::PyExpr::new(op.density_matrix().clone()))
        );

        impl_expr_getter_interface!(
            fn_name = $overlap_distribution_fn,
            fn_doc = impl_expr_getter_doc!("overlap distribution", $expr_ty),
            expr_ty = $expr_ty,
            out_ty = $crate::core::expr::PyExpr,
            body =
                |op: &$expr_ty| Ok($crate::core::expr::PyExpr::new(op.overlap_distribution().clone()))
        );

        impl_expr_getter_interface!(
            fn_name = $grid_expr_fn,
            fn_doc = impl_expr_getter_doc!("grid expression", $expr_ty),
            expr_ty = $expr_ty,
            out_ty = $crate::core::expr::PyExpr,
            body = |op: &$expr_ty| Ok($crate::core::expr::PyExpr::new(op.$grid_expr_method().clone()))
        );

        impl_expr_getter_interface!(
            fn_name = $derivative_fn,
            fn_doc = impl_expr_getter_doc!("derivative", $expr_ty),
            expr_ty = $expr_ty,
            out_ty = $crate::perturbations::pert_multichain::PyPertMultichain,
            body = |op: &$expr_ty| Ok($crate::perturbations::pert_multichain::PyPertMultichain::new(
                op.derivative().clone(),
            ))
        );

        pub fn $register_fn(
            m: &::pyo3::Bound<'_, ::pyo3::types::PyModule>,
        ) -> ::pyo3::PyResult<()> {
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
