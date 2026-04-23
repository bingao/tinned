use pyo3::prelude::*;

use tinned::AoTwoElecMatrix;

use crate::core::{errors::to_pyerr, expr::PyExpr};
use crate::perturbations::pert_multichain::PyPertMultichain;

/// Create an AoTwoElecMatrix expression.
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
#[pyo3(signature = (name, density, dependencies=None, derivative=None))]
pub fn ao_two_elec_matrix_new(
    name: String,
    density: &PyExpr,
    dependencies: Option<&PyPertMultichain>,
    derivative: Option<&PyPertMultichain>,
) -> PyResult<PyExpr> {
    let mut b = AoTwoElecMatrix::builder(name, density.inner().clone());

    if let Some(v) = dependencies {
        b = b.dependencies(v.inner().clone());
    }
    if let Some(v) = derivative {
        b = b.derivative(v.inner().clone());
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = ao_two_elec_matrix_name,
    fn_doc = impl_expr_getter_doc!("name", AoTwoElecMatrix),
    expr_ty = AoTwoElecMatrix,
    out_ty = String,
    body = |op: &AoTwoElecMatrix| Ok(op.name().to_string())
);

impl_expr_getter_interface!(
    fn_name = ao_two_elec_matrix_density,
    fn_doc = impl_expr_getter_doc!("density", AoTwoElecMatrix),
    expr_ty = AoTwoElecMatrix,
    out_ty = PyExpr,
    body = |op: &AoTwoElecMatrix| Ok(PyExpr::new(op.density().clone()))
);

impl_expr_getter_interface!(
    fn_name = ao_two_elec_matrix_dependencies,
    fn_doc = impl_expr_getter_doc!("dependencies", AoTwoElecMatrix),
    expr_ty = AoTwoElecMatrix,
    out_ty = PyPertMultichain,
    body = |op: &AoTwoElecMatrix| Ok(PyPertMultichain::new(op.dependencies().clone()))
);

impl_expr_getter_interface!(
    fn_name = ao_two_elec_matrix_derivative,
    fn_doc = impl_expr_getter_doc!("derivative", AoTwoElecMatrix),
    expr_ty = AoTwoElecMatrix,
    out_ty = PyPertMultichain,
    body = |op: &AoTwoElecMatrix| Ok(PyPertMultichain::new(op.derivative().clone()))
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ao_two_elec_matrix_new, m)?)?;
    m.add_function(wrap_pyfunction!(ao_two_elec_matrix_name, m)?)?;
    m.add_function(wrap_pyfunction!(ao_two_elec_matrix_density, m)?)?;
    m.add_function(wrap_pyfunction!(ao_two_elec_matrix_dependencies, m)?)?;
    m.add_function(wrap_pyfunction!(ao_two_elec_matrix_derivative, m)?)?;
    Ok(())
}
