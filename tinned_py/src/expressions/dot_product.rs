use pyo3::prelude::*;
use std::sync::Arc;

use tinned::{DotProduct, Expr, TinnedError};

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
pub fn dot_product_new(
    bra: PyExpr,
    use_hermitian: bool,
    ket: PyExpr,
    allow_braket_swap: bool,
    is_scalar: Option<bool>,
) -> PyResult<PyExpr> {
    let rust_bra: Arc<dyn Expr> = bra.inner().clone();
    let rust_ket: Arc<dyn Expr> = ket.inner().clone();

    let out = DotProduct::new(rust_bra, use_hermitian, rust_ket, allow_braket_swap, is_scalar)
        .map_err(to_pyerr)?;

    Ok(PyExpr::new(out))
}

/// Return the bra of a DotProduct expression.
#[pyfunction]
pub fn dot_product_bra(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let dp_ref = inner.as_any().downcast_ref::<DotProduct>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "dot_product_bra() expected a DotProduct expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(dp_ref.bra().clone()))
}

/// Return the ket of a DotProduct expression.
#[pyfunction]
pub fn dot_product_ket(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let dp_ref = inner.as_any().downcast_ref::<DotProduct>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "dot_product_ket() expected a DotProduct expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(dp_ref.ket().clone()))
}

/// Return whether bra-ket swapping is allowed for this DotProduct.
#[pyfunction]
pub fn dot_product_allow_braket_swap(expr: PyExpr) -> PyResult<bool> {
    let inner = expr.inner().clone();

    let dp_ref = inner.as_any().downcast_ref::<DotProduct>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "dot_product_allow_braket_swap() expected a DotProduct expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(dp_ref.allow_braket_swap())
}

/// Return the conjugate of a DotProduct expression.
#[pyfunction]
pub fn dot_product_conjugate(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let dp_ref = inner.as_any().downcast_ref::<DotProduct>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "dot_product_conjugate() expected a DotProduct expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    let result = dp_ref.conjugate().map_err(to_pyerr)?;
    Ok(PyExpr::new(result))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dot_product_new, m)?)?;
    m.add_function(wrap_pyfunction!(dot_product_bra, m)?)?;
    m.add_function(wrap_pyfunction!(dot_product_ket, m)?)?;
    m.add_function(wrap_pyfunction!(dot_product_allow_braket_swap, m)?)?;
    m.add_function(wrap_pyfunction!(dot_product_conjugate, m)?)?;
    Ok(())
}
