use pyo3::prelude::*;

use tinned::{TinnedError, TwoElecOperator};

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

/// Return name for a TwoElecOperator.
///
/// Errors if the input expression is not a TwoElecOperator.
#[pyfunction]
pub fn two_elec_operator_name(expr: PyExpr) -> PyResult<String> {
    let inner = expr.inner().clone();

    let op_ref = inner.as_any().downcast_ref::<TwoElecOperator>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "two_elec_operator_name() expected a TwoElecOperator expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(op_ref.name().to_string())
}

/// Return density for a TwoElecOperator as a PyExpr.
///
/// Errors if the input expression is not a TwoElecOperator.
#[pyfunction]
pub fn two_elec_operator_density(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let op_ref = inner.as_any().downcast_ref::<TwoElecOperator>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "two_elec_operator_density() expected a TwoElecOperator expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(op_ref.density().clone()))
}

/// Return dependencies for a TwoElecOperator.
///
/// Errors if the input expression is not a TwoElecOperator.
#[pyfunction]
pub fn two_elec_operator_dependencies(expr: PyExpr) -> PyResult<PyPertMultichain> {
    let inner = expr.inner().clone();

    let op_ref = inner.as_any().downcast_ref::<TwoElecOperator>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "two_elec_operator_dependencies() expected a TwoElecOperator expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyPertMultichain::new(op_ref.dependencies().clone()))
}

/// Return derivative for a TwoElecOperator.
///
/// Errors if the input expression is not a TwoElecOperator.
#[pyfunction]
pub fn two_elec_operator_derivative(expr: PyExpr) -> PyResult<PyPertMultichain> {
    let inner = expr.inner().clone();

    let op_ref = inner.as_any().downcast_ref::<TwoElecOperator>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "two_elec_operator_derivative() expected a TwoElecOperator expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyPertMultichain::new(op_ref.derivative().clone()))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(two_elec_operator_new, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_operator_name, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_operator_density, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_operator_dependencies, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_operator_derivative, m)?)?;
    Ok(())
}
