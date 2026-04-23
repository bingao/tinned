macro_rules! impl_unary_expr_interface {
    (
        expr_ty = $expr_ty:ty,
        new_fn = $new_fn:ident,
        new_doc = $new_doc:literal,
        argument_fn = $argument_fn:ident,
        argument_doc = $argument_doc:literal,
        register_fn = $register_fn:ident
    ) => {
        #[doc = $new_doc]
        #[::pyo3::prelude::pyfunction]
        pub fn $new_fn(argument: py_expr_ref_ty!()) -> ::pyo3::PyResult<py_expr_ty!()> {
            let rust_argument: ::std::sync::Arc<dyn ::tinned::Expr> = argument.inner().clone();
            let out = <$expr_ty>::new(rust_argument).map_err($crate::core::errors::to_pyerr)?;
            Ok(<py_expr_ty!()>::new(out))
        }

        impl_expr_getter_interface!(
            fn_name = $argument_fn,
            fn_doc = $argument_doc,
            expr_ty = $expr_ty,
            out_ty = py_expr_ty!(),
            body = |op: &$expr_ty| Ok(<py_expr_ty!()>::new(op.argument().clone()))
        );

        pub fn $register_fn(
            m: &::pyo3::Bound<'_, ::pyo3::types::PyModule>,
        ) -> ::pyo3::PyResult<()> {
            m.add_function(::pyo3::wrap_pyfunction!($new_fn, m)?)?;
            m.add_function(::pyo3::wrap_pyfunction!($argument_fn, m)?)?;
            Ok(())
        }
    };
}
