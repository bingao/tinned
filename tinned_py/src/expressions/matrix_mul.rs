use pyo3::prelude::*;
use std::vec::Vec;

use tinned::MatrixMul;

use crate::core::{errors::to_pyerr, expr::PyExpr};

/// Build a matrix multiplication from terms.
///
/// Args:
///   terms: Terms (scalars become part of the coefficient; non-scalars are factors).
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
pub fn matrix_mul_new(terms: Vec<PyExpr>) -> PyResult<PyExpr> {
    let rust_terms = terms.into_iter().map(|t| t.inner().clone()).collect();

    let out = MatrixMul::new(rust_terms).map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = matrix_mul_coefficient,
    fn_doc = impl_expr_getter_doc!("coefficient", MatrixMul),
    expr_ty = MatrixMul,
    out_ty = PyExpr,
    body = |mul: &MatrixMul| Ok(PyExpr::new(mul.coefficient().clone()))
);

impl_expr_getter_interface!(
    fn_name = matrix_mul_factors,
    fn_doc = impl_expr_getter_doc!("factors (excluding the coefficient)", MatrixMul),
    expr_ty = MatrixMul,
    out_ty = Vec<PyExpr>,
    body = |mul: &MatrixMul| Ok(mul.factors().iter().cloned().map(PyExpr::new).collect())
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(matrix_mul_new, m)?)?;
    m.add_function(wrap_pyfunction!(matrix_mul_coefficient, m)?)?;
    m.add_function(wrap_pyfunction!(matrix_mul_factors, m)?)?;

    Ok(())
}
