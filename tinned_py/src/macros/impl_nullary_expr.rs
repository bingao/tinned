macro_rules! impl_nullary_expr_interface {
    (
        expr_ty = $expr_ty:ty,
        has_deps = $has_deps:tt,
        new_fn = $new_fn:ident,
        name_fn = $name_fn:ident,
        is_perturbing_fn = $is_perturbing_fn:ident,
        derivative_fn = $derivative_fn:ident,
        $(
            deps_fn = $deps_fn:ident,
            indep_perts_fn = $indep_perts_fn:ident,
        )?
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
            fn_name = $is_perturbing_fn,
            fn_doc = impl_expr_getter_doc!("Whether the operator is a perturbing one or not", $expr_ty),
            expr_ty = $expr_ty,
            out_ty = bool,
            body = |op: &$expr_ty| Ok(op.is_perturbing())
        );

        impl_expr_getter_interface!(
            fn_name = $derivative_fn,
            fn_doc = impl_expr_getter_doc!("derivative", $expr_ty),
            expr_ty = $expr_ty,
            out_ty = py_pert_multichain_ty!(),
            body = |op: &$expr_ty| Ok(<py_pert_multichain_ty!()>::new(
                op.derivative().clone(),
            ))
        );

        $(
            impl_nullary_expr_interface!(@with_deps_fn
                expr_ty = $expr_ty,
                deps_fn = $deps_fn,
                indep_perts_fn = $indep_perts_fn
            );
        )?

        pub fn $register_fn(m: &::pyo3::Bound<'_, ::pyo3::types::PyModule>) -> ::pyo3::PyResult<()> {
            m.add_function(::pyo3::wrap_pyfunction!($new_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($name_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($is_perturbing_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($derivative_fn, m)?)?;
            $(
                m.add_function(::pyo3::wrap_pyfunction!($deps_fn, m)?)?;
                m.add_function(::pyo3::wrap_pyfunction!($indep_perts_fn, m)?)?;
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
            "  is_perturbing: Optional bool indicating a perturbing expression or not (default False).\n",
            "  dependencies: Optional PertMultichain (default empty).\n",
            "  independent_perturbations: Optional no mixed derivatives are allowed within independent perturbations (default empty).\n",
            "  derivative: Optional PertMultichain (default empty).\n\n",
            "Returns:\n",
            "  A PyExpr wrapping the constructed expression (interned).\n\n",
            "Notes:\n",
            "  If dependencies is not a superchain of derivative, the Rust constructor returns a zero expression."
        )]
        #[::pyo3::prelude::pyfunction]
        #[pyo3(signature = (
            name,
            is_perturbing=None,
            dependencies=None,
            independent_perturbations=None,
            derivative=None
        ))]
        pub fn $new_fn(
            name: ::std::string::String,
            is_perturbing: ::std::option::Option<bool>,
            dependencies: ::std::option::Option<py_pert_multichain_ref_ty!()>,
            independent_perturbations: ::std::option::Option<::std::vec::Vec<py_pert_ty!()>>,
            derivative: ::std::option::Option<py_pert_multichain_ref_ty!()>,
        ) -> ::pyo3::PyResult<py_expr_ty!()> {
            let mut b = <$expr_ty>::builder(name);

            if let Some(v) = is_perturbing {
                b = b.is_perturbing(v);
            }
            if let Some(deps) = dependencies {
                b = b.dependencies(deps.inner().clone());
            }
            if let Some(indep_perts) = independent_perturbations {
                let indep_perts = indep_perts
                    .into_iter()
                    .map(|p| p.inner().clone())
                    .collect::<::std::collections::BTreeSet<_>>();

                b = b.independent_perturbations(indep_perts);
            }
            if let Some(deriv) = derivative {
                b = b.derivative(deriv.inner().clone());
            }

            let out = b.build().map_err($crate::core::errors::to_pyerr)?;
            Ok(<py_expr_ty!()>::new(out))
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
            "  is_perturbing: Optional bool indicating a perturbing expression or not (default False).\n",
            "  derivative: Optional PertMultichain (default empty).\n\n",
            "Returns:\n",
            "  A PyExpr wrapping the constructed expression (interned)."
        )]
        #[::pyo3::prelude::pyfunction]
        #[pyo3(signature = (
            name,
            is_perturbing=None,
            derivative=None
        ))]
        pub fn $new_fn(
            name: ::std::string::String,
            is_perturbing: ::std::option::Option<bool>,
            derivative: ::std::option::Option<py_pert_multichain_ref_ty!()>,
        ) -> ::pyo3::PyResult<py_expr_ty!()> {
            let mut b = <$expr_ty>::builder(name);

            if let Some(v) = is_perturbing {
                b = b.is_perturbing(v);
            }
            if let Some(deriv) = derivative {
                b = b.derivative(deriv.inner().clone());
            }

            let out = b.build().map_err($crate::core::errors::to_pyerr)?;
            Ok(<py_expr_ty!()>::new(out))
        }
    };

    (@with_deps_fn
        expr_ty = $expr_ty:ty,
        deps_fn = $deps_fn:ident,
        indep_perts_fn = $indep_perts_fn:ident
    ) => {
        impl_expr_getter_interface!(
            fn_name = $deps_fn,
            fn_doc = impl_expr_getter_doc!("dependencies", $expr_ty),
            expr_ty = $expr_ty,
            out_ty = py_pert_multichain_ty!(),
            body = |op: &$expr_ty| Ok(<py_pert_multichain_ty!()>::new(
                op.dependencies().clone(),
            ))
        );

        impl_expr_getter_interface!(
            fn_name = $indep_perts_fn,
            fn_doc = impl_expr_getter_doc!("independent perturbations", $expr_ty),
            expr_ty = $expr_ty,
            out_ty = ::std::vec::Vec<py_pert_ty!()>,
            body = |op: &$expr_ty| Ok(
                op.independent_perturbations()
                    .iter()
                    .cloned()
                    .map(<py_pert_ty!()>::new)
                    .collect::<Vec<_>>()
            )
        );
    };
}
