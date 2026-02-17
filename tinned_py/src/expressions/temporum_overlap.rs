use pyo3::prelude::*;

use tinned::{TemporumOverlap, TinnedError};

use crate::core::{errors::to_pyerr, expr::PyExpr};
use crate::perturbations::pert_multichain::PyPertMultichain;

/// Create a TemporumOverlap expression.
///
/// Args:
///   dependencies: Perturbation dependencies of Sb and Sk.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
pub fn temporum_overlap_new(dependencies: &Bound<'_, PyPertMultichain>) -> PyResult<PyExpr> {
    let deps = dependencies.borrow().inner().clone();

    let b = TemporumOverlap::builder(deps);
    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Return is_zero_strength for a TemporumOverlap.
///
/// Errors if the input expression is not a TemporumOverlap.
#[pyfunction]
pub fn temporum_overlap_is_zero_strength(expr: PyExpr) -> PyResult<bool> {
    let inner = expr.inner().clone();

    let ov_ref = inner.as_any().downcast_ref::<TemporumOverlap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "temporum_overlap_is_zero_strength() expected a TemporumOverlap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(ov_ref.is_zero_strength())
}

/// Return braket for a TemporumOverlap as a PyExpr.
///
/// Errors if the input expression is not a TemporumOverlap.
#[pyfunction]
pub fn temporum_overlap_braket(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let ov_ref = inner.as_any().downcast_ref::<TemporumOverlap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "temporum_overlap_braket() expected a TemporumOverlap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(ov_ref.braket().clone()))
}

/// Return dependencies for a TemporumOverlap.
///
/// Errors if the input expression is not a TemporumOverlap.
#[pyfunction]
pub fn temporum_overlap_dependencies(expr: PyExpr) -> PyResult<PyPertMultichain> {
    let inner = expr.inner().clone();

    let ov_ref = inner.as_any().downcast_ref::<TemporumOverlap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "temporum_overlap_dependencies() expected a TemporumOverlap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyPertMultichain::new(ov_ref.dependencies().clone()))
}

/// Return derivative for a TemporumOverlap.
///
/// Errors if the input expression is not a TemporumOverlap.
#[pyfunction]
pub fn temporum_overlap_derivative(expr: PyExpr) -> PyResult<PyPertMultichain> {
    let inner = expr.inner().clone();

    let ov_ref = inner.as_any().downcast_ref::<TemporumOverlap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "temporum_overlap_derivative() expected a TemporumOverlap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyPertMultichain::new(ov_ref.derivative().clone()))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(temporum_overlap_new, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_overlap_is_zero_strength, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_overlap_braket, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_overlap_dependencies, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_overlap_derivative, m)?)?;
    Ok(())
}
