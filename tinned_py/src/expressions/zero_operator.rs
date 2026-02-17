use pyo3::prelude::*;
use std::sync::Arc;

use tinned::{Expr, ZeroOperator};

use crate::core::expr::PyExpr;

/// Build a ZeroOperator.
///
/// Returns:
///   A PyExpr wrapping the zero operator (interned).
#[pyfunction]
pub fn zero_operator_new() -> PyResult<PyExpr> {
    let out: Arc<dyn Expr> = ZeroOperator::new();
    Ok(PyExpr::new(out))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(zero_operator_new, m)?)?;
    Ok(())
}
