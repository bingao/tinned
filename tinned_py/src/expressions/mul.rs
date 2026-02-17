use pyo3::prelude::*;
use std::sync::Arc;
use std::vec::Vec;

use tinned::{Expr, Mul, TinnedError};

use crate::core::{errors::to_pyerr, expr::PyExpr};

/// Build a multiplication from scalar terms.
///
/// Args:
///   terms: Scalar terms.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
pub fn mul_new(terms: Vec<PyExpr>) -> PyResult<PyExpr> {
    let rust_terms: Vec<Arc<dyn Expr>> = terms.into_iter().map(|t| t.inner().clone()).collect();

    let out = Mul::new(rust_terms).map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Return the coefficient of a Mul as a PyExpr.
///
/// Errors if the input expression is not a Mul.
#[pyfunction]
pub fn mul_coefficient(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let mul_ref = inner.as_any().downcast_ref::<Mul>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "mul_coefficient() expected a Mul expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    // Convert &Number to Arc<dyn Expr> and wrap in PyExpr.
    let coeff_expr: Arc<dyn Expr> = mul_ref.coefficient().clone().into();

    Ok(PyExpr::new(coeff_expr))
}

/// Return the factors of a Mul (excluding the numeric coefficient).
///
/// Errors if the input expression is not a Mul.
#[pyfunction]
pub fn mul_factors(expr: PyExpr) -> PyResult<Vec<PyExpr>> {
    let inner = expr.inner().clone();

    let mul_ref = inner.as_any().downcast_ref::<Mul>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "mul_factors() expected a Mul expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(mul_ref.factors().iter().cloned().map(PyExpr::new).collect())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(mul_new, m)?)?;
    m.add_function(wrap_pyfunction!(mul_coefficient, m)?)?;
    m.add_function(wrap_pyfunction!(mul_factors, m)?)?;

    Ok(())
}
