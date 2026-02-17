macro_rules! impl_addition_interface {
    (
        new_fn = $new_fn:ident,
        new_fn_doc = $new_fn_doc:literal,
        terms_fn = $terms_fn:ident,
        terms_fn_doc = $terms_fn_doc:literal,
        register_fn = $register_fn:ident,
        expr_ty = $expr_ty:ty,
        downcast_err_msg = $downcast_err_msg:literal
    ) => {
        #[doc = $new_fn_doc]
        #[pyo3::prelude::pyfunction]
        pub fn $new_fn(
            terms: ::std::vec::Vec<$crate::core::expr::PyExpr>,
        ) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let rust_terms: ::std::vec::Vec<::std::sync::Arc<dyn tinned::Expr>> =
                terms.into_iter().map(|t| t.inner().clone()).collect();

            let out = <$expr_ty>::new(rust_terms).map_err($crate::core::errors::to_pyerr)?;
            Ok($crate::core::expr::PyExpr::new(out))
        }

        #[doc = $terms_fn_doc]
        #[pyo3::prelude::pyfunction]
        pub fn $terms_fn(
            expr: $crate::core::expr::PyExpr,
        ) -> ::pyo3::PyResult<::std::vec::Vec<$crate::core::expr::PyExpr>> {
            let inner = expr.inner().clone();

            let add_ref = inner.as_any().downcast_ref::<$expr_ty>().ok_or_else(|| {
                $crate::core::errors::to_pyerr(tinned::TinnedError::ExpressionError {
                    message: $downcast_err_msg,
                    expression: inner.to_string(),
                    source: None,
                })
            })?;

            Ok(add_ref.terms().iter().cloned().map($crate::core::expr::PyExpr::new).collect())
        }

        pub fn $register_fn(
            m: &::pyo3::Bound<'_, ::pyo3::types::PyModule>,
        ) -> ::pyo3::PyResult<()> {
            m.add_function(::pyo3::wrap_pyfunction!($new_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($terms_fn, m)?)?;
            Ok(())
        }
    };
}
