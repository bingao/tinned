use pyo3::prelude::*;
use std::sync::Arc;

use tinned::{Composition, Expr};

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

impl_expr_getter_interface!(
    fn_name = composition_name,
    fn_doc = impl_expr_getter_doc!("name", Composition),
    expr_ty = Composition,
    out_ty = String,
    body = |op: &Composition| Ok(op.name().to_string())
);

impl_expr_getter_interface!(
    fn_name = composition_order,
    fn_doc = impl_expr_getter_doc!("order", Composition),
    expr_ty = Composition,
    out_ty = u32,
    body = |op: &Composition| Ok(op.order())
);

impl_expr_getter_interface!(
    fn_name = composition_inner,
    fn_doc = impl_expr_getter_doc!("inner expression", Composition),
    expr_ty = Composition,
    out_ty = PyExpr,
    body = |op: &Composition| Ok(PyExpr::new(op.inner().clone()))
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(composition_new, m)?)?;
    m.add_function(wrap_pyfunction!(composition_name, m)?)?;
    m.add_function(wrap_pyfunction!(composition_order, m)?)?;
    m.add_function(wrap_pyfunction!(composition_inner, m)?)?;
    Ok(())
}
