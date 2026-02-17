use pyo3::prelude::*;

use tinned::{ExpAdjointMap, TinnedError};

use crate::core::{errors::to_pyerr, expr::PyExpr};
use crate::perturbations::pert_multichain::PyPertMultichain;

/// Create an ExpAdjointMap expression.
///
/// Args:
///   generator: Generator expression.
///   target: Target expression.
///   left_action: Optional bool.
///                If True: exp(X)*Y*exp(-X).
///                If False: exp(-X)*Y*exp(X).
///   max_fold: Optional u32 truncation.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
pub fn exp_adjoint_map_new(
    generator: PyExpr,
    target: PyExpr,
    left_action: Option<bool>,
    max_fold: Option<u32>,
) -> PyResult<PyExpr> {
    let mut b = ExpAdjointMap::builder(generator.inner().clone(), target.inner().clone());

    if let Some(v) = left_action {
        b = b.left_action(v);
    }
    if let Some(v) = max_fold {
        b = b.max_fold(v);
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Create an ExpAdjointMap expression in temporum mode.
///
/// Args:
///   generator: Generator expression.
///   is_forward: If True uses i*d/dt, otherwise -i*d/dt.
///   left_action: Optional bool.
///   max_fold: Optional u32 truncation.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
pub fn exp_adjoint_map_temporum_new(
    generator: PyExpr,
    is_forward: bool,
    left_action: Option<bool>,
    max_fold: Option<u32>,
) -> PyResult<PyExpr> {
    let mut b = ExpAdjointMap::builder_temporum(generator.inner().clone(), is_forward);

    if let Some(v) = left_action {
        b = b.left_action(v);
    }
    if let Some(v) = max_fold {
        b = b.max_fold(v);
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Return the generator of an ExpAdjointMap expression.
#[pyfunction]
pub fn exp_adjoint_map_generator(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let m_ref = inner.as_any().downcast_ref::<ExpAdjointMap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "exp_adjoint_map_generator() expected an ExpAdjointMap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(m_ref.generator().clone()))
}

/// Return the target of an ExpAdjointMap expression.
#[pyfunction]
pub fn exp_adjoint_map_target(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let m_ref = inner.as_any().downcast_ref::<ExpAdjointMap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "exp_adjoint_map_target() expected an ExpAdjointMap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(m_ref.target().clone()))
}

/// Return whether this ExpAdjointMap is in temporum mode.
///
/// If True, target is obtained by performing time differentiation
/// on generator.
#[pyfunction]
pub fn exp_adjoint_map_is_temporum(expr: PyExpr) -> PyResult<bool> {
    let inner = expr.inner().clone();

    let m_ref = inner.as_any().downcast_ref::<ExpAdjointMap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "exp_adjoint_map_is_temporum() expected an ExpAdjointMap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(m_ref.is_temporum())
}

/// Return whether this ExpAdjointMap uses left action.
///
/// If True: exp(X)*Y*exp(-X).
/// If False: exp(-X)*Y*exp(X).
#[pyfunction]
pub fn exp_adjoint_map_left_action(expr: PyExpr) -> PyResult<bool> {
    let inner = expr.inner().clone();

    let m_ref = inner.as_any().downcast_ref::<ExpAdjointMap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "exp_adjoint_map_left_action() expected an ExpAdjointMap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(m_ref.left_action())
}

/// Return the maximum fold of an ExpAdjointMap expression.
#[pyfunction]
pub fn exp_adjoint_map_max_fold(expr: PyExpr) -> PyResult<u32> {
    let inner = expr.inner().clone();

    let m_ref = inner.as_any().downcast_ref::<ExpAdjointMap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "exp_adjoint_map_max_fold() expected an ExpAdjointMap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(m_ref.max_fold())
}

/// Return whether the generator strength is zero for this ExpAdjointMap expression.
#[pyfunction]
pub fn exp_adjoint_map_is_zero_strength(expr: PyExpr) -> PyResult<bool> {
    let inner = expr.inner().clone();

    let m_ref = inner.as_any().downcast_ref::<ExpAdjointMap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "exp_adjoint_map_is_zero_strength() expected an ExpAdjointMap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(m_ref.is_zero_strength())
}

/// Return the result expression stored inside ExpAdjointMap.
///
/// This is the differentiated expression of the exponential adjoint map.
#[pyfunction]
pub fn exp_adjoint_map_result(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let m_ref = inner.as_any().downcast_ref::<ExpAdjointMap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "exp_adjoint_map_result() expected an ExpAdjointMap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(m_ref.result().clone()))
}

/// Return the derivative of an ExpAdjointMap expression.
#[pyfunction]
pub fn exp_adjoint_map_derivative(expr: PyExpr) -> PyResult<PyPertMultichain> {
    let inner = expr.inner().clone();

    let m_ref = inner.as_any().downcast_ref::<ExpAdjointMap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "exp_adjoint_map_derivative() expected an ExpAdjointMap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyPertMultichain::new(m_ref.derivative().clone()))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(exp_adjoint_map_new, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_temporum_new, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_generator, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_target, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_is_temporum, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_left_action, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_max_fold, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_is_zero_strength, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_result, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_derivative, m)?)?;

    Ok(())
}
