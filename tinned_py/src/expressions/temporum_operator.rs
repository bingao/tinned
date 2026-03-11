use pyo3::prelude::*;

use tinned::{TemporumOperator, TinnedError, generic_expression_error};

use crate::core::{errors::to_pyerr, expr::PyExpr};
use crate::perturbations::pert_multichain::PyPertMultichain;

/// Create a TemporumOperator expression.
///
/// Args:
///   argument: Non-scalar argument expression.
///   is_forward: Optional bool. If True uses i*d/dt, otherwise -i*d/dt.
///               Defaults to True.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
pub fn temporum_operator_new(argument: PyExpr, is_forward: Option<bool>) -> PyResult<PyExpr> {
    let mut b = TemporumOperator::builder(argument.inner().clone());

    if let Some(v) = is_forward {
        b = b.is_forward(v);
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Return is_forward for a TemporumOperator.
///
/// Errors if the input expression is not a TemporumOperator.
#[pyfunction]
pub fn temporum_operator_is_forward(expr: PyExpr) -> PyResult<bool> {
    let inner = expr.inner().clone();

    let top_ref = inner.as_any().downcast_ref::<TemporumOperator>().ok_or_else(|| {
to_pyerr(generic_expression_error(
    "temporum_operator_is_forward() expected a TemporumOperator expression",
    &inner,
    None,
))    })?;

    Ok(top_ref.is_forward())
}

/// Return argument for a TemporumOperator as a PyExpr.
///
/// Errors if the input expression is not a TemporumOperator.
#[pyfunction]
pub fn temporum_operator_argument(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let top_ref = inner.as_any().downcast_ref::<TemporumOperator>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "temporum_operator_argument() expected a TemporumOperator expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(top_ref.argument().clone()))
}

/// Return derivative for a TemporumOperator.
///
/// Errors if the input expression is not a TemporumOperator.
#[pyfunction]
pub fn temporum_operator_derivative(expr: PyExpr) -> PyResult<PyPertMultichain> {
    let inner = expr.inner().clone();

    let top_ref = inner.as_any().downcast_ref::<TemporumOperator>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "temporum_operator_derivative() expected a TemporumOperator expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    let d = top_ref.derivative().map_err(to_pyerr)?;
    Ok(PyPertMultichain::new(d.clone()))
}

/// Return frequency expression for a TemporumOperator.
///
/// For unperturbed argument, frequency returns zero.
///
/// Errors if the input expression is not a TemporumOperator.
#[pyfunction]
pub fn temporum_operator_frequency(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let top_ref = inner.as_any().downcast_ref::<TemporumOperator>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "temporum_operator_frequency() expected a TemporumOperator expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    let f = top_ref.frequency().map_err(to_pyerr)?;
    Ok(PyExpr::new(f))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(temporum_operator_new, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_operator_is_forward, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_operator_argument, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_operator_derivative, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_operator_frequency, m)?)?;

    Ok(())
}
