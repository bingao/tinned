macro_rules! impl_addition_interface {
    (
        expr_ty = $expr_ty:ty,
        new_fn = $new_fn:ident,
        new_doc = $new_doc:literal,
        terms_fn = $terms_fn:ident,
        terms_doc = $terms_doc:literal,
        register_fn = $register_fn:ident
    ) => {
        #[doc = $new_doc]
        #[::pyo3::prelude::pyfunction]
        pub fn $new_fn(
            terms: ::std::vec::Vec<py_expr_ty!()>,
        ) -> ::pyo3::PyResult<py_expr_ty!()> {
            let rust_terms: ::std::vec::Vec<::std::sync::Arc<dyn ::tinned::Expr>> =
                terms.into_iter().map(|t| t.inner().clone()).collect();

            let out = <$expr_ty>::new(rust_terms).map_err($crate::core::errors::to_pyerr)?;
            Ok(<py_expr_ty!()>::new(out))
        }

        impl_expr_getter_interface!(
            fn_name = $terms_fn,
            fn_doc = $terms_doc,
            expr_ty = $expr_ty,
            out_ty = ::std::vec::Vec<py_expr_ty!()>,
            body = |add: &$expr_ty| Ok(add.terms().iter().cloned().map(<py_expr_ty!()>::new).collect())
        );

        pub fn $register_fn(
            m: &::pyo3::Bound<'_, ::pyo3::types::PyModule>,
        ) -> ::pyo3::PyResult<()> {
            m.add_function(::pyo3::wrap_pyfunction!($new_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($terms_fn, m)?)?;
            Ok(())
        }
    };
}
