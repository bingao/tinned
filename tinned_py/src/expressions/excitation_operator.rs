use pyo3::prelude::*;
use std::sync::Arc;

use tinned::{ExcitationOperator, Expr};

use crate::core::expr::PyExpr;

/// Build an excitation operator with a given name.
///
/// Args:
///   name: name of the excitation operator.
///
/// Returns:
///   A PyExpr wrapping the constructed excitation operator (interned).
#[pyfunction]
pub fn excitation_operator_new(name: String) -> PyResult<PyExpr> {
    let out: Arc<dyn Expr> = ExcitationOperator::new(name);
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = excitation_operator_name,
    fn_doc = impl_expr_getter_doc!("name", ExcitationOperator),
    expr_ty = ExcitationOperator,
    out_ty = String,
    body = |op: &ExcitationOperator| Ok(op.name().to_string())
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(excitation_operator_new, m)?)?;
    m.add_function(wrap_pyfunction!(excitation_operator_name, m)?)?;

    Ok(())
}
