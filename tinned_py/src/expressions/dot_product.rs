use pyo3::prelude::*;

use tinned::DotProduct;

use crate::core::{errors::to_pyerr, expr::PyExpr};

/// Create a DotProduct expression.
///
/// Args:
///   bra: Bra expression (non-scalar).
///   use_hermitian: If True, apply HermitianTranspose to bra;
///                  otherwise Transpose.
///   ket: Ket expression (non-scalar).
///   allow_braket_swap: If True, allow canonical swapping of bra and ket.
///   is_scalar: Optional scalar flag (default True).
///
/// Returns:
///   A PyExpr wrapping the constructed DotProduct (interned).
///
/// Notes:
///   - bra and ket must be non-scalar.
///   - If either side is ZeroOperator, the result is zero.
#[pyfunction]
#[pyo3(signature = (
    bra,
    use_hermitian,
    ket,
    allow_braket_swap,
    is_scalar=true
))]
pub fn dot_product_new(
    bra: &PyExpr,
    use_hermitian: bool,
    ket: &PyExpr,
    allow_braket_swap: bool,
    is_scalar: bool,
) -> PyResult<PyExpr> {
    let rust_bra = bra.inner().clone();
    let rust_ket = ket.inner().clone();

    let out =
        DotProduct::new(rust_bra, use_hermitian, rust_ket, allow_braket_swap, Some(is_scalar))
            .map_err(to_pyerr)?;

    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = dot_product_bra,
    fn_doc = impl_expr_getter_doc!("bra", DotProduct),
    expr_ty = DotProduct,
    out_ty = PyExpr,
    body = |op: &DotProduct| Ok(PyExpr::new(op.bra().clone()))
);

impl_expr_getter_interface!(
    fn_name = dot_product_ket,
    fn_doc = impl_expr_getter_doc!("ket", DotProduct),
    expr_ty = DotProduct,
    out_ty = PyExpr,
    body = |op: &DotProduct| Ok(PyExpr::new(op.ket().clone()))
);

impl_expr_getter_interface!(
    fn_name = dot_product_allow_braket_swap,
    fn_doc = impl_expr_getter_doc!("whether bra-ket swapping is allowed", DotProduct),
    expr_ty = DotProduct,
    out_ty = bool,
    body = |op: &DotProduct| Ok(op.allow_braket_swap())
);

impl_expr_getter_interface!(
    fn_name = dot_product_conjugate,
    fn_doc = impl_expr_getter_doc!("conjugate", DotProduct),
    expr_ty = DotProduct,
    out_ty = PyExpr,
    body = |op: &DotProduct| {
        let result = op.conjugate().map_err(to_pyerr)?;
        Ok(PyExpr::new(result))
    }
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dot_product_new, m)?)?;
    m.add_function(wrap_pyfunction!(dot_product_bra, m)?)?;
    m.add_function(wrap_pyfunction!(dot_product_ket, m)?)?;
    m.add_function(wrap_pyfunction!(dot_product_allow_braket_swap, m)?)?;
    m.add_function(wrap_pyfunction!(dot_product_conjugate, m)?)?;
    Ok(())
}
