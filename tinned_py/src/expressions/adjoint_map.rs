use pyo3::prelude::*;
use std::sync::Arc;

use tinned::{AdjointMap, Expr};

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
#[pyo3(signature = (generators, target, left_action=true))]
pub fn adjoint_map_new(
    generators: Vec<PyExpr>,
    target: PyExpr,
    left_action: bool,
) -> PyResult<PyExpr> {
    let rust_generators: Vec<Arc<dyn Expr>> =
        generators.into_iter().map(|t| t.inner().clone()).collect();
    let rust_target: Arc<dyn Expr> = target.inner().clone();

    let out = AdjointMap::new(rust_generators, rust_target, Some(left_action)).map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = adjoint_map_generators,
    fn_doc = impl_expr_getter_doc!("generators", AdjointMap),
    expr_ty = AdjointMap,
    out_ty = Vec<PyExpr>,
    body = |op: &AdjointMap| Ok(op.generators().iter().cloned().map(PyExpr::new).collect())
);

impl_expr_getter_interface!(
    fn_name = adjoint_map_target,
    fn_doc = impl_expr_getter_doc!("target", AdjointMap),
    expr_ty = AdjointMap,
    out_ty = PyExpr,
    body = |op: &AdjointMap| Ok(PyExpr::new(op.target().clone()))
);

impl_expr_getter_interface!(
    fn_name = adjoint_map_left_action,
    fn_doc = impl_expr_getter_doc!("Boolean value of left action", AdjointMap),
    expr_ty = AdjointMap,
    out_ty = bool,
    body = |op: &AdjointMap| Ok(op.left_action())
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(adjoint_map_new, m)?)?;
    m.add_function(wrap_pyfunction!(adjoint_map_generators, m)?)?;
    m.add_function(wrap_pyfunction!(adjoint_map_target, m)?)?;
    m.add_function(wrap_pyfunction!(adjoint_map_left_action, m)?)?;
    Ok(())
}
