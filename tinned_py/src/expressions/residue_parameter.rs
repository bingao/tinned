use pyo3::prelude::*;
use std::vec::Vec;

use tinned::ResidueParameter;

use crate::core::{errors::to_pyerr, expr::PyExpr};
use crate::perturbations::{pert_multichain::PyPertMultichain, perturbation::PyPerturbation};

/// Create a ResidueParameter expression.
///
/// Args:
///   perturbations: Iterable of Perturbation.
///   excited_state: Excited state expression.
///   parameter: Perturbed parameter expression.
///   positive_frequency: Optional bool which indicates the sum of frequencies of
///      some perturbations approaches the energy of an excited state from the
///      positive or negative side. Defaults to True, positive side.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned), or ZeroOperator under the builder rules.
#[pyfunction]
#[pyo3(signature = (perturbations, excited_state, parameter, positive_frequency=None))]
pub fn residue_parameter_new(
    perturbations: Vec<PyPerturbation>,
    excited_state: PyExpr,
    parameter: PyExpr,
    positive_frequency: Option<bool>,
) -> PyResult<PyExpr> {
    let perts = perturbations.into_iter().map(|p| p.inner().clone()).collect();

    let mut b =
        ResidueParameter::builder(perts, excited_state.inner().clone(), parameter.inner().clone());

    if let Some(v) = positive_frequency {
        b = b.positive_frequency(v);
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = residue_parameter_positive_frequency,
    fn_doc = impl_expr_getter_doc!(
        "whether the sum of perturbation frequencies approaches an exitation energy from the positive side",
        ResidueParameter
    ),
    expr_ty = ResidueParameter,
    out_ty = bool,
    body = |op: &ResidueParameter| Ok(op.positive_frequency())
);

impl_expr_getter_interface!(
    fn_name = residue_parameter_perturbations,
    fn_doc = impl_expr_getter_doc!("perturbations", ResidueParameter),
    expr_ty = ResidueParameter,
    out_ty = Vec<PyPerturbation>,
    body = |residue: &ResidueParameter| Ok(residue.perturbations().iter().cloned().map(PyPerturbation::new).collect())
);

impl_expr_getter_interface!(
    fn_name = residue_parameter_excited_state,
    fn_doc = impl_expr_getter_doc!("excited state", ResidueParameter),
    expr_ty = ResidueParameter,
    out_ty = PyExpr,
    body = |residue: &ResidueParameter| Ok(PyExpr::new(residue.excited_state().clone()))
);

impl_expr_getter_interface!(
    fn_name = residue_parameter_parameter,
    fn_doc = impl_expr_getter_doc!("perturbed parameter", ResidueParameter),
    expr_ty = ResidueParameter,
    out_ty = PyExpr,
    body = |residue: &ResidueParameter| Ok(PyExpr::new(residue.parameter().clone()))
);

impl_expr_getter_interface!(
    fn_name = residue_parameter_derivative,
    fn_doc = impl_expr_getter_doc!("derivative", ResidueParameter),
    expr_ty = ResidueParameter,
    out_ty = PyPertMultichain,
    body = |residue: &ResidueParameter| {
        let derivative = residue.derivative().map_err(to_pyerr)?;
        Ok(PyPertMultichain::new(derivative.clone()))
    }
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(residue_parameter_new, m)?)?;
    m.add_function(wrap_pyfunction!(residue_parameter_positive_frequency, m)?)?;
    m.add_function(wrap_pyfunction!(residue_parameter_perturbations, m)?)?;
    m.add_function(wrap_pyfunction!(residue_parameter_excited_state, m)?)?;
    m.add_function(wrap_pyfunction!(residue_parameter_parameter, m)?)?;
    m.add_function(wrap_pyfunction!(residue_parameter_derivative, m)?)?;

    Ok(())
}
