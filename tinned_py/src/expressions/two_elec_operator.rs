use pyo3::prelude::*;

use tinned::TwoElecOperator;

use crate::core::{errors::to_pyerr, expr::PyExpr};
use crate::perturbations::pert_multichain::PyPertMultichain;

/// Create a TwoElecOperator expression.
///
/// Args:
///   name: Operator name.
///   density: Density expression (must be WfnParameter or ResidueParameter; ZeroOperator passes through).
///   dependencies: Optional dependencies. Defaults to empty.
///   derivative: Optional derivative. Defaults to empty.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned), or ZeroOperator when derivative is not a subchain of dependencies.
#[pyfunction]
pub fn two_elec_operator_new(
    name: String,
    density: PyExpr,
    dependencies: Option<&Bound<'_, PyPertMultichain>>,
    derivative: Option<&Bound<'_, PyPertMultichain>>,
) -> PyResult<PyExpr> {
    let mut b = TwoElecOperator::builder(name, density.inner().clone());

    if let Some(v) = dependencies {
        b = b.dependencies(v.borrow().inner().clone());
    }
    if let Some(v) = derivative {
        b = b.derivative(v.borrow().inner().clone());
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = two_elec_operator_name,
    fn_doc = impl_expr_getter_doc!("name", TwoElecOperator),
    expr_ty = TwoElecOperator,
    out_ty = String,
    body = |op: &TwoElecOperator| Ok(op.name().to_string())
);

impl_expr_getter_interface!(
    fn_name = two_elec_operator_density,
    fn_doc = impl_expr_getter_doc!("density", TwoElecOperator),
    expr_ty = TwoElecOperator,
    out_ty = PyExpr,
    body = |op: &TwoElecOperator| Ok(PyExpr::new(op.density().clone()))
);

impl_expr_getter_interface!(
    fn_name = two_elec_operator_dependencies,
    fn_doc = impl_expr_getter_doc!("dependencies", TwoElecOperator),
    expr_ty = TwoElecOperator,
    out_ty = PyPertMultichain,
    body = |op: &TwoElecOperator| Ok(PyPertMultichain::new(op.dependencies().clone()))
);

impl_expr_getter_interface!(
    fn_name = two_elec_operator_derivative,
    fn_doc = impl_expr_getter_doc!("derivative", TwoElecOperator),
    expr_ty = TwoElecOperator,
    out_ty = PyPertMultichain,
    body = |op: &TwoElecOperator| Ok(PyPertMultichain::new(op.derivative().clone()))
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(two_elec_operator_new, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_operator_name, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_operator_density, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_operator_dependencies, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_operator_derivative, m)?)?;
    Ok(())
}
