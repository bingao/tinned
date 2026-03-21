use pyo3::prelude::*;

use tinned::TemporumOverlap;

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
pub fn temporum_overlap_new(dependencies: &PyPertMultichain) -> PyResult<PyExpr> {
    let deps = dependencies.inner().clone();

    let b = TemporumOverlap::builder(deps);
    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = temporum_overlap_zero_rules_applied,
    fn_doc = impl_expr_getter_doc!("whether evaluated at zero-field strength", TemporumOverlap),
    expr_ty = TemporumOverlap,
    out_ty = bool,
    body = |op: &TemporumOverlap| Ok(op.zero_rules_applied())
);

impl_expr_getter_interface!(
    fn_name = temporum_overlap_braket,
    fn_doc = impl_expr_getter_doc!("braket expression", TemporumOverlap),
    expr_ty = TemporumOverlap,
    out_ty = PyExpr,
    body = |op: &TemporumOverlap| Ok(PyExpr::new(op.braket().clone()))
);

impl_expr_getter_interface!(
    fn_name = temporum_overlap_dependencies,
    fn_doc = impl_expr_getter_doc!("dependencies", TemporumOverlap),
    expr_ty = TemporumOverlap,
    out_ty = PyPertMultichain,
    body = |op: &TemporumOverlap| Ok(PyPertMultichain::new(op.dependencies().clone()))
);

impl_expr_getter_interface!(
    fn_name = temporum_overlap_derivative,
    fn_doc = impl_expr_getter_doc!("derivative", TemporumOverlap),
    expr_ty = TemporumOverlap,
    out_ty = PyPertMultichain,
    body = |op: &TemporumOverlap| Ok(PyPertMultichain::new(op.derivative().clone()))
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(temporum_overlap_new, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_overlap_zero_rules_applied, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_overlap_braket, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_overlap_dependencies, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_overlap_derivative, m)?)?;
    Ok(())
}
