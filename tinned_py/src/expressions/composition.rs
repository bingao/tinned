use pyo3::prelude::*;
use std::sync::Arc;

use tinned::{Composition, Expr, TinnedError};

use crate::core::{errors::to_pyerr, expr::PyExpr};

/// Create a Composition expression.
///
/// Args:
///   name: Name of the outer function.
///   order: Order parameter.
///   inner: A scalar expression.
///
/// Returns:
///   A PyExpr wrapping the constructed Composition (interned).
///
/// Notes:
///   - inner must be scalar.
///   - If inner is zero, the result is Number::zero().
#[pyfunction]
pub fn composition_new(name: String, order: u32, inner: PyExpr) -> PyResult<PyExpr> {
    let rust_inner: Arc<dyn Expr> = inner.inner().clone();

    let out = Composition::new(name, order, rust_inner).map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Return the name of a Composition expression.
#[pyfunction]
pub fn composition_name(expr: PyExpr) -> PyResult<String> {
    let inner = expr.inner().clone();

    let comp_ref = inner.as_any().downcast_ref::<Composition>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "composition_name() expected a Composition expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(comp_ref.name().to_string())
}

/// Return the order of a Composition expression.
#[pyfunction]
pub fn composition_order(expr: PyExpr) -> PyResult<u32> {
    let inner = expr.inner().clone();

    let comp_ref = inner.as_any().downcast_ref::<Composition>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "composition_order() expected a Composition expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(comp_ref.order())
}

/// Return the inner expression of a Composition expression.
#[pyfunction]
pub fn composition_inner(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let comp_ref = inner.as_any().downcast_ref::<Composition>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "composition_inner() expected a Composition expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(comp_ref.inner().clone()))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(composition_new, m)?)?;
    m.add_function(wrap_pyfunction!(composition_name, m)?)?;
    m.add_function(wrap_pyfunction!(composition_order, m)?)?;
    m.add_function(wrap_pyfunction!(composition_inner, m)?)?;
    Ok(())
}
