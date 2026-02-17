macro_rules! impl_nullary_expr_interface {
    (
        type_name = $type_name:ty,
        type_label = $type_label:literal,
        has_deps = $has_deps:tt,
        new_fn = $new_fn:ident,
        name_fn = $name_fn:ident,
        derivative_fn = $derivative_fn:ident,
        $( deps_fn = $deps_fn:ident, )?
        register_fn = $register_fn:ident
    ) => {
        impl_nullary_expr_interface!(@new_fn
            type_name = $type_name,
            type_label = $type_label,
            has_deps = $has_deps,
            new_fn = $new_fn
        );

        #[doc = concat!("Return the name of a ", $type_label, " expression.")]
        #[pyo3::prelude::pyfunction]
        pub fn $name_fn(expr: $crate::core::expr::PyExpr) -> ::pyo3::PyResult<::std::string::String> {
            let inner = expr.inner().clone();

            let op_ref = inner.as_any().downcast_ref::<$type_name>().ok_or_else(|| {
                $crate::core::errors::to_pyerr(tinned::TinnedError::ExpressionError {
                    message: concat!(stringify!($name_fn), "() expected a ", $type_label, " expression"),
                    expression: inner.to_string(),
                    source: None,
                })
            })?;

            Ok(op_ref.name().to_string())
        }

        #[doc = concat!("Return the derivative of a ", $type_label, " expression.")]
        #[pyo3::prelude::pyfunction]
        pub fn $derivative_fn(
            expr: $crate::core::expr::PyExpr,
        ) -> ::pyo3::PyResult<$crate::perturbations::pert_multichain::PyPertMultichain> {
            let inner = expr.inner().clone();

            let op_ref = inner.as_any().downcast_ref::<$type_name>().ok_or_else(|| {
                $crate::core::errors::to_pyerr(tinned::TinnedError::ExpressionError {
                    message: concat!(stringify!($derivative_fn), "() expected a ", $type_label, " expression"),
                    expression: inner.to_string(),
                    source: None,
                })
            })?;

            Ok($crate::perturbations::pert_multichain::PyPertMultichain::new(
                op_ref.derivative().clone(),
            ))
        }

        $(
            impl_nullary_expr_interface!(@deps_fn
                type_name = $type_name,
                type_label = $type_label,
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
        type_name = $type_name:ty,
        type_label = $type_label:literal,
        has_deps = true,
        new_fn = $new_fn:ident
    ) => {
        #[doc = concat!(
            "Create a ", $type_label, " expression.\n\n",
            "Args:\n",
            "  name: Name.\n",
            "  derivative: Optional PertMultichain (default empty).\n",
            "  dependencies: Optional PertMultichain (default empty).\n\n",
            "Returns:\n",
            "  A PyExpr wrapping the constructed expression (interned).\n\n",
            "Notes:\n",
            "  If dependencies is not a superchain of derivative, the Rust constructor returns a zero expression."
        )]
        #[pyo3::prelude::pyfunction]
        pub fn $new_fn(
            name: ::std::string::String,
            derivative: ::std::option::Option<&::pyo3::Bound<'_, $crate::perturbations::pert_multichain::PyPertMultichain>>,
            dependencies: ::std::option::Option<&::pyo3::Bound<'_, $crate::perturbations::pert_multichain::PyPertMultichain>>,
        ) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let mut b = <$type_name>::builder(name);

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
        type_name = $type_name:ty,
        type_label = $type_label:literal,
        has_deps = false,
        new_fn = $new_fn:ident
    ) => {
        #[doc = concat!(
            "Create a ", $type_label, " expression.\n\n",
            "Args:\n",
            "  name: Name.\n",
            "  derivative: Optional PertMultichain (default empty).\n\n",
            "Returns:\n",
            "  A PyExpr wrapping the constructed expression (interned)."
        )]
        #[pyo3::prelude::pyfunction]
        pub fn $new_fn(
            name: ::std::string::String,
            derivative: ::std::option::Option<&::pyo3::Bound<'_, $crate::perturbations::pert_multichain::PyPertMultichain>>,
        ) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let mut b = <$type_name>::builder(name);

            if let Some(deriv) = derivative {
                b = b.derivative(deriv.borrow().inner().clone());
            }

            let out = b.build().map_err($crate::core::errors::to_pyerr)?;
            Ok($crate::core::expr::PyExpr::new(out))
        }
    };

    (@deps_fn
        type_name = $type_name:ty,
        type_label = $type_label:literal,
        deps_fn = $deps_fn:ident
    ) => {
        #[doc = concat!("Return the dependencies of a ", $type_label, " expression.")]
        #[pyo3::prelude::pyfunction]
        pub fn $deps_fn(
            expr: $crate::core::expr::PyExpr,
        ) -> ::pyo3::PyResult<$crate::perturbations::pert_multichain::PyPertMultichain> {
            let inner = expr.inner().clone();

            let op_ref = inner.as_any().downcast_ref::<$type_name>().ok_or_else(|| {
                $crate::core::errors::to_pyerr(tinned::TinnedError::ExpressionError {
                    message: concat!(stringify!($deps_fn), "() expected a ", $type_label, " expression"),
                    expression: inner.to_string(),
                    source: None,
                })
            })?;

            Ok($crate::perturbations::pert_multichain::PyPertMultichain::new(
                op_ref.dependencies().clone(),
            ))
        }
    };
}
