use pyo3::prelude::*;
use std::sync::Arc;

use tinned::{Expr, Power, TinnedError};

use crate::core::{errors::to_pyerr, expr::PyExpr};

/// Build a power expression from base and exponent.
///
/// Args:
///   base: Scalar base expression.
///   exponent: Integer exponent.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
pub fn power_new(base: PyExpr, exponent: i64) -> PyResult<PyExpr> {
    let base_inner: Arc<dyn Expr> = base.inner().clone();

    let out = Power::new(base_inner, exponent).map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Return the base of a Power as a PyExpr.
///
/// Errors if the input expression is not a Power.
#[pyfunction]
pub fn power_base(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let pow_ref = inner.as_any().downcast_ref::<Power>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "power_base() expected a Power expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(pow_ref.base().clone()))
}

/// Return the exponent of a Power.
///
/// Errors if the input expression is not a Power.
#[pyfunction]
pub fn power_exponent(expr: PyExpr) -> PyResult<i64> {
    let inner = expr.inner().clone();

    let pow_ref = inner.as_any().downcast_ref::<Power>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "power_exponent() expected a Power expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(pow_ref.exponent())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(power_new, m)?)?;
    m.add_function(wrap_pyfunction!(power_base, m)?)?;
    m.add_function(wrap_pyfunction!(power_exponent, m)?)?;

    Ok(())
}
