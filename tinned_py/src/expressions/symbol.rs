use pyo3::prelude::*;
use std::sync::Arc;

use tinned::{Expr, Symbol, TinnedError};

use crate::core::{errors::to_pyerr, expr::PyExpr};

/// Build a Symbol from a name.
///
/// Args:
///   name: Symbol name.
///
/// Returns:
///   A PyExpr wrapping the constructed symbol (interned).
#[pyfunction]
pub fn symbol_new(name: String) -> PyResult<PyExpr> {
    let out: Arc<dyn Expr> = Symbol::new(name);
    Ok(PyExpr::new(out))
}

/// Return the name of a Symbol.
///
/// Errors if the input expression is not a Symbol.
#[pyfunction]
pub fn symbol_name(expr: PyExpr) -> PyResult<String> {
    let inner = expr.inner().clone();

    let sym_ref = inner.as_any().downcast_ref::<Symbol>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "symbol_name() expected a Symbol expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(sym_ref.name().to_string())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(symbol_new, m)?)?;
    m.add_function(wrap_pyfunction!(symbol_name, m)?)?;

    Ok(())
}
