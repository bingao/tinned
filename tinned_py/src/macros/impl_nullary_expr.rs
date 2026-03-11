macro_rules! impl_nullary_expr_interface {
    (
        expr_ty = $expr_ty:ty,
        has_deps = $has_deps:tt,
        new_fn = $new_fn:ident,
        name_fn = $name_fn:ident,
        derivative_fn = $derivative_fn:ident,
        $( deps_fn = $deps_fn:ident, )?
        register_fn = $register_fn:ident
    ) => {
        impl_nullary_expr_interface!(@new_fn
            expr_ty = $expr_ty,
            has_deps = $has_deps,
            new_fn = $new_fn
        );

        impl_expr_getter_interface!(
            fn_name = $name_fn,
            fn_doc = impl_expr_getter_doc!("name", $expr_ty),
            expr_ty = $expr_ty,
            out_ty = ::std::string::String,
            body = |op: &$expr_ty| Ok(op.name().to_string())
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

        $(
            impl_nullary_expr_interface!(@deps_fn
                expr_ty = $expr_ty,
                deps_fn = $deps_fn
            );
        )?

        pub fn $register_fn(m: &::pyo3::Bound<'_, ::pyo3::types::PyModule>) -> ::pyo3::PyResult<()> {
            m.add_function(::pyo3::wrap_pyfunction!($new_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($name_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($derivative_fn, m)?)?;
            $(
                m.add_function(::pyo3::wrap_pyfunction!($deps_fn, m)?)?;
            )?
            Ok(())
        }
    };

    (@new_fn
        expr_ty = $expr_ty:ty,
        has_deps = true,
        new_fn = $new_fn:ident
    ) => {
        #[doc = concat!(
            "Create a ", stringify!($expr_ty), " expression.\n\n",
            "Args:\n",
            "  name: Name.\n",
            "  derivative: Optional PertMultichain (default empty).\n",
            "  dependencies: Optional PertMultichain (default empty).\n\n",
            "Returns:\n",
            "  A PyExpr wrapping the constructed expression (interned).\n\n",
            "Notes:\n",
            "  If dependencies is not a superchain of derivative, the Rust constructor returns a zero expression."
        )]
        #[::pyo3::prelude::pyfunction]
        pub fn $new_fn(
            name: ::std::string::String,
            derivative: ::std::option::Option<&::pyo3::Bound<'_, $crate::perturbations::pert_multichain::PyPertMultichain>>,
            dependencies: ::std::option::Option<&::pyo3::Bound<'_, $crate::perturbations::pert_multichain::PyPertMultichain>>,
        ) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let mut b = <$expr_ty>::builder(name);

            if let Some(deps) = dependencies {
                b = b.dependencies(deps.borrow().inner().clone());
            }
            if let Some(deriv) = derivative {
                b = b.derivative(deriv.borrow().inner().clone());
            }

            let out = b.build().map_err($crate::core::errors::to_pyerr)?;
            Ok($crate::core::expr::PyExpr::new(out))
        }
    };

    (@new_fn
        expr_ty = $expr_ty:ty,
        has_deps = false,
        new_fn = $new_fn:ident
    ) => {
        #[doc = concat!(
            "Create a ", stringify!($expr_ty), " expression.\n\n",
            "Args:\n",
            "  name: Name.\n",
            "  derivative: Optional PertMultichain (default empty).\n\n",
            "Returns:\n",
            "  A PyExpr wrapping the constructed expression (interned)."
        )]
        #[::pyo3::prelude::pyfunction]
        pub fn $new_fn(
            name: ::std::string::String,
            derivative: ::std::option::Option<&::pyo3::Bound<'_, $crate::perturbations::pert_multichain::PyPertMultichain>>,
        ) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let mut b = <$expr_ty>::builder(name);

            if let Some(deriv) = derivative {
                b = b.derivative(deriv.borrow().inner().clone());
            }

            let out = b.build().map_err($crate::core::errors::to_pyerr)?;
            Ok($crate::core::expr::PyExpr::new(out))
        }
    };

    (@deps_fn
        expr_ty = $expr_ty:ty,
        deps_fn = $deps_fn:ident
    ) => {
        impl_expr_getter_interface!(
            fn_name = $deps_fn,
            fn_doc = impl_expr_getter_doc!("dependencies", $expr_ty),
            expr_ty = $expr_ty,
            out_ty = $crate::perturbations::pert_multichain::PyPertMultichain,
            body = |op: &$expr_ty| Ok($crate::perturbations::pert_multichain::PyPertMultichain::new(
                op.dependencies().clone(),
            ))
        );
    };
}
