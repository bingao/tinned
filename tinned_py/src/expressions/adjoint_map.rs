use pyo3::prelude::*;
use std::sync::Arc;

use tinned::{AdjointMap, Expr, TinnedError};

use crate::core::{errors::to_pyerr, expr::PyExpr};

/// Create an AdjointMap expression.
///
/// Args:
///   generators: A list of non-scalar operator expressions.
///   target: A non-scalar operator expression.
///   left_action: If provided, sets whether the adjoint map acts
///                from the left (default True).
///
/// Returns:
///   A PyExpr wrapping the constructed AdjointMap (interned).
#[pyfunction]
pub fn adjoint_map_new(
    generators: Vec<PyExpr>,
    target: PyExpr,
    left_action: Option<bool>,
) -> PyResult<PyExpr> {
    let rust_generators: Vec<Arc<dyn Expr>> =
        generators.into_iter().map(|t| t.inner().clone()).collect();
    let rust_target: Arc<dyn Expr> = target.inner().clone();

    let out = AdjointMap::new(rust_generators, rust_target, left_action).map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Return the generators of an AdjointMap expression as a list.
#[pyfunction]
pub fn adjoint_map_generators(expr: PyExpr) -> PyResult<Vec<PyExpr>> {
    let inner = expr.inner().clone();

    let adj_ref = inner.as_any().downcast_ref::<AdjointMap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "adjoint_map_generators() expected an AdjointMap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(adj_ref.generators().iter().cloned().map(PyExpr::new).collect())
}

/// Return the target of an AdjointMap expression.
#[pyfunction]
pub fn adjoint_map_target(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let adj_ref = inner.as_any().downcast_ref::<AdjointMap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "adjoint_map_target() expected an AdjointMap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(adj_ref.target().clone()))
}

/// Return whether an AdjointMap expression uses left action.
#[pyfunction]
pub fn adjoint_map_left_action(expr: PyExpr) -> PyResult<bool> {
    let inner = expr.inner().clone();

    let adj_ref = inner.as_any().downcast_ref::<AdjointMap>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "adjoint_map_left_action() expected an AdjointMap expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(adj_ref.left_action())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(adjoint_map_new, m)?)?;
    m.add_function(wrap_pyfunction!(adjoint_map_generators, m)?)?;
    m.add_function(wrap_pyfunction!(adjoint_map_target, m)?)?;
    m.add_function(wrap_pyfunction!(adjoint_map_left_action, m)?)?;
    Ok(())
}
