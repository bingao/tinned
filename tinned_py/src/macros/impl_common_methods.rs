macro_rules! impl_expr_getter_interface {
    (
        fn_name = $fn_name:ident,
        fn_doc = $fn_doc:expr,
        expr_ty = $expr_ty:ty,
        out_ty = $out_ty:ty,
        body = $body:expr
    ) => {
        #[doc = $fn_doc]
        #[::pyo3::prelude::pyfunction]
        pub fn $fn_name(expr: py_expr_ref_ty!()) -> ::pyo3::PyResult<$out_ty> {
            let inner = expr.inner().clone();

            let op_ref = inner.as_any().downcast_ref::<$expr_ty>().ok_or_else(|| {
                $crate::core::errors::to_pyerr(::tinned::public::expression_error(
                    concat!(
                        stringify!($fn_name),
                        "() expected a ",
                        stringify!($expr_ty),
                        " expression"
                    ),
                    &inner,
                    None,
                ))
            })?;

            ($body)(op_ref)
        }
    };
}

macro_rules! impl_expr_getter_doc {
    ($field:literal, $type:ty) => {
        concat!(
            "Return ",
            $field,
            " for an expression ",
            stringify!($type),
            ".\n\n",
            "Raises an error if the expression is not ",
            stringify!($type),
            "."
        )
    };
}
