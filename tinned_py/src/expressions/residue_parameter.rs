use pyo3::prelude::*;
use std::sync::Arc;
use std::vec::Vec;

use tinned::{ResidueParameter, TinnedError};

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
///      positive or negative side. Defaults to True.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned), or ZeroOperator under the builder rules.
#[pyfunction]
pub fn residue_parameter_new(
    perturbations: &Bound<'_, PyAny>,
    excited_state: PyExpr,
    parameter: PyExpr,
    positive_frequency: Option<bool>,
) -> PyResult<PyExpr> {
    let mut perts: Vec<Arc<tinned::Perturbation>> = Vec::new();

    for item in perturbations.try_iter()? {
        let item = item?;
        let p: Bound<'_, PyPerturbation> = item.extract()?;
        perts.push(p.borrow().inner().clone());
    }

    let mut b =
        ResidueParameter::builder(perts, excited_state.inner().clone(), parameter.inner().clone());

    if let Some(v) = positive_frequency {
        b = b.positive_frequency(v);
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Return whether the sum of frequencies of some perturbations approaches the
/// energy of an excited state from the positive or negative side.
///
/// Errors if the input expression is not a ResidueParameter.
#[pyfunction]
pub fn residue_parameter_positive_frequency(expr: PyExpr) -> PyResult<bool> {
    let inner = expr.inner().clone();

    let res_ref = inner
        .as_any()
        .downcast_ref::<ResidueParameter>()
        .ok_or_else(|| {
            to_pyerr(TinnedError::ExpressionError {
                message: "residue_parameter_positive_frequency() expected a ResidueParameter expression",
                expression: inner.to_string(),
                source: None,
            })
        })?;

    Ok(res_ref.positive_frequency())
}

/// Return perturbations for a ResidueParameter.
///
/// Errors if the input expression is not a ResidueParameter.
#[pyfunction]
pub fn residue_parameter_perturbations(expr: PyExpr) -> PyResult<Vec<PyPerturbation>> {
    let inner = expr.inner().clone();

    let res_ref = inner.as_any().downcast_ref::<ResidueParameter>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "residue_parameter_perturbations() expected a ResidueParameter expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(res_ref.perturbations().iter().cloned().map(PyPerturbation::new).collect())
}

/// Return excited state for a ResidueParameter as a PyExpr.
///
/// Errors if the input expression is not a ResidueParameter.
#[pyfunction]
pub fn residue_parameter_excited_state(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let res_ref = inner.as_any().downcast_ref::<ResidueParameter>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "residue_parameter_excited_state() expected a ResidueParameter expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(res_ref.excited_state().clone()))
}

/// Return perturbed parameter for a ResidueParameter as a PyExpr.
///
/// Errors if the input expression is not a ResidueParameter.
#[pyfunction]
pub fn residue_parameter_parameter(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let res_ref = inner.as_any().downcast_ref::<ResidueParameter>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "residue_parameter_parameter() expected a ResidueParameter expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(res_ref.parameter().clone()))
}

/// Return derivative for a ResidueParameter.
///
/// Errors if the input expression is not a ResidueParameter.
#[pyfunction]
pub fn residue_parameter_derivative(expr: PyExpr) -> PyResult<PyPertMultichain> {
    let inner = expr.inner().clone();

    let res_ref = inner.as_any().downcast_ref::<ResidueParameter>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "residue_parameter_derivative() expected a ResidueParameter expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    let d = res_ref.derivative().map_err(to_pyerr)?;
    Ok(PyPertMultichain::new(d.clone()))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(residue_parameter_new, m)?)?;
    m.add_function(wrap_pyfunction!(residue_parameter_positive_frequency, m)?)?;
    m.add_function(wrap_pyfunction!(residue_parameter_perturbations, m)?)?;
    m.add_function(wrap_pyfunction!(residue_parameter_excited_state, m)?)?;
    m.add_function(wrap_pyfunction!(residue_parameter_parameter, m)?)?;
    m.add_function(wrap_pyfunction!(residue_parameter_derivative, m)?)?;

    Ok(())
}
