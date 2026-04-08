use pyo3::prelude::*;

use tinned::BasisTimeEvolution;

use crate::core::{errors::to_pyerr, expr::PyExpr};
use crate::perturbations::pert_multichain::PyPertMultichain;

/// Create a BasisTimeEvolution expression.
///
/// Args:
///   dependencies: Perturbation dependencies of Sb and Sk.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
pub fn basis_time_evolution_new(dependencies: &PyPertMultichain) -> PyResult<PyExpr> {
    let deps = dependencies.inner().clone();

    let b = BasisTimeEvolution::builder(deps);
    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = basis_time_evolution_at_zero_perturbations,
    fn_doc = impl_expr_getter_doc!("whether evaluated at zero-field strength", BasisTimeEvolution),
    expr_ty = BasisTimeEvolution,
    out_ty = bool,
    body = |op: &BasisTimeEvolution| Ok(op.at_zero_perturbations())
);

impl_expr_getter_interface!(
    fn_name = basis_time_evolution_braket,
    fn_doc = impl_expr_getter_doc!("braket expression", BasisTimeEvolution),
    expr_ty = BasisTimeEvolution,
    out_ty = PyExpr,
    body = |op: &BasisTimeEvolution| Ok(PyExpr::new(op.braket().clone()))
);

impl_expr_getter_interface!(
    fn_name = basis_time_evolution_dependencies,
    fn_doc = impl_expr_getter_doc!("dependencies", BasisTimeEvolution),
    expr_ty = BasisTimeEvolution,
    out_ty = PyPertMultichain,
    body = |op: &BasisTimeEvolution| Ok(PyPertMultichain::new(op.dependencies().clone()))
);

impl_expr_getter_interface!(
    fn_name = basis_time_evolution_derivative,
    fn_doc = impl_expr_getter_doc!("derivative", BasisTimeEvolution),
    expr_ty = BasisTimeEvolution,
    out_ty = PyPertMultichain,
    body = |op: &BasisTimeEvolution| Ok(PyPertMultichain::new(op.derivative().clone()))
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(basis_time_evolution_new, m)?)?;
    m.add_function(wrap_pyfunction!(basis_time_evolution_at_zero_perturbations, m)?)?;
    m.add_function(wrap_pyfunction!(basis_time_evolution_braket, m)?)?;
    m.add_function(wrap_pyfunction!(basis_time_evolution_dependencies, m)?)?;
    m.add_function(wrap_pyfunction!(basis_time_evolution_derivative, m)?)?;
    Ok(())
}
