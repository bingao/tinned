macro_rules! impl_unary_expr_interface {
    (
        new_fn = $new_fn:ident,
        new_fn_doc = $new_fn_doc:literal,
        arg_fn = $arg_fn:ident,
        arg_fn_doc = $arg_fn_doc:literal,
        register_fn = $register_fn:ident,
        expr_ty = $expr_ty:ty,
        downcast_err_msg = $downcast_err_msg:literal
    ) => {
        #[doc = $new_fn_doc]
        #[pyo3::prelude::pyfunction]
        pub fn $new_fn(
            argument: $crate::core::expr::PyExpr,
        ) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let rust_argument: ::std::sync::Arc<dyn tinned::Expr> = argument.inner().clone();
            let out = <$expr_ty>::new(rust_argument).map_err($crate::core::errors::to_pyerr)?;
            Ok($crate::core::expr::PyExpr::new(out))
        }

        #[doc = $arg_fn_doc]
        #[pyo3::prelude::pyfunction]
        pub fn $arg_fn(
            expr: $crate::core::expr::PyExpr,
        ) -> ::pyo3::PyResult<$crate::core::expr::PyExpr> {
            let inner = expr.inner().clone();

            let op_ref = inner.as_any().downcast_ref::<$expr_ty>().ok_or_else(|| {
                $crate::core::errors::to_pyerr(tinned::TinnedError::ExpressionError {
                    message: $downcast_err_msg,
                    expression: inner.to_string(),
                    source: None,
                })
            })?;

            Ok($crate::core::expr::PyExpr::new(op_ref.argument().clone()))
        }

        pub fn $register_fn(
            m: &::pyo3::Bound<'_, ::pyo3::types::PyModule>,
        ) -> ::pyo3::PyResult<()> {
            m.add_function(::pyo3::wrap_pyfunction!($new_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($arg_fn, m)?)?;
            Ok(())
        }
    };
}
