use pyo3::prelude::*;

use tinned::{TinnedError, TwoElecEnergy, TwoElecOperator};

use crate::core::{errors::to_pyerr, expr::PyExpr};
use crate::perturbations::pert_multichain::PyPertMultichain;

/// Create a TwoElecEnergy expression.
///
/// Args:
///   name: Energy name.
///   inner_density: Inner density expression (must be WfnParameter or ResidueParameter; ZeroOperator yields Number(0)).
///   outer_density: Optional outer density expression. Defaults to inner_density.
///   allow_density_swap: Optional bool controlling equality behavior. Defaults to True.
///   dependencies: Optional dependencies. Defaults to empty.
///   derivative: Optional derivative. Defaults to empty.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned), or Number(0) when invalid/zero per rules.
#[pyfunction]
pub fn two_elec_energy_new(
    name: String,
    inner_density: PyExpr,
    outer_density: Option<PyExpr>,
    allow_density_swap: Option<bool>,
    dependencies: Option<&Bound<'_, PyPertMultichain>>,
    derivative: Option<&Bound<'_, PyPertMultichain>>,
) -> PyResult<PyExpr> {
    let mut b = TwoElecEnergy::builder(name, inner_density.inner().clone());

    if let Some(v) = outer_density {
        b = b.outer_density(v.inner().clone());
    }
    if let Some(v) = allow_density_swap {
        b = b.allow_density_swap(v);
    }
    if let Some(v) = dependencies {
        b = b.dependencies(v.borrow().inner().clone());
    }
    if let Some(v) = derivative {
        b = b.derivative(v.borrow().inner().clone());
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Create a TwoElecEnergy expression from a TwoElecOperator expression.
///
/// Args:
///   two_elec_op: A TwoElecOperator expression.
///   outer_density: Optional outer density expression. Defaults to operator density.
///   allow_density_swap: Optional bool controlling equality behavior. Defaults to True.
///   dependencies: Optional dependencies. Defaults to operator dependencies.
///   derivative: Optional derivative. Defaults to operator derivative.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned), or Number(0) when invalid/zero per rules.
#[pyfunction]
pub fn two_elec_energy_from_operator(
    two_elec_op: PyExpr,
    outer_density: Option<PyExpr>,
    allow_density_swap: Option<bool>,
    dependencies: Option<&Bound<'_, PyPertMultichain>>,
    derivative: Option<&Bound<'_, PyPertMultichain>>,
) -> PyResult<PyExpr> {
    let inner = two_elec_op.inner().clone();

    let op_ref = inner.as_any().downcast_ref::<TwoElecOperator>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "two_elec_energy_from_operator() expected a TwoElecOperator expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    let mut b = TwoElecEnergy::builder_from_operator(op_ref);

    if let Some(v) = outer_density {
        b = b.outer_density(v.inner().clone());
    }
    if let Some(v) = allow_density_swap {
        b = b.allow_density_swap(v);
    }
    if let Some(v) = dependencies {
        b = b.dependencies(v.borrow().inner().clone());
    }
    if let Some(v) = derivative {
        b = b.derivative(v.borrow().inner().clone());
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Return name for a TwoElecEnergy.
///
/// Errors if the input expression is not a TwoElecEnergy.
#[pyfunction]
pub fn two_elec_energy_name(expr: PyExpr) -> PyResult<String> {
    let inner = expr.inner().clone();

    let e_ref = inner.as_any().downcast_ref::<TwoElecEnergy>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "two_elec_energy_name() expected a TwoElecEnergy expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(e_ref.name().to_string())
}

/// Return inner_density for a TwoElecEnergy as a PyExpr.
///
/// Errors if the input expression is not a TwoElecEnergy.
#[pyfunction]
pub fn two_elec_energy_inner_density(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let e_ref = inner.as_any().downcast_ref::<TwoElecEnergy>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "two_elec_energy_inner_density() expected a TwoElecEnergy expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(e_ref.inner_density().clone()))
}

/// Return outer_density for a TwoElecEnergy as a PyExpr.
///
/// Errors if the input expression is not a TwoElecEnergy.
#[pyfunction]
pub fn two_elec_energy_outer_density(expr: PyExpr) -> PyResult<PyExpr> {
    let inner = expr.inner().clone();

    let e_ref = inner.as_any().downcast_ref::<TwoElecEnergy>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "two_elec_energy_outer_density() expected a TwoElecEnergy expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyExpr::new(e_ref.outer_density().clone()))
}

/// Return allow_density_swap for a TwoElecEnergy.
///
/// Errors if the input expression is not a TwoElecEnergy.
#[pyfunction]
pub fn two_elec_energy_allow_density_swap(expr: PyExpr) -> PyResult<bool> {
    let inner = expr.inner().clone();

    let e_ref = inner.as_any().downcast_ref::<TwoElecEnergy>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "two_elec_energy_allow_density_swap() expected a TwoElecEnergy expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(e_ref.allow_density_swap())
}

/// Return dependencies for a TwoElecEnergy.
///
/// Errors if the input expression is not a TwoElecEnergy.
#[pyfunction]
pub fn two_elec_energy_dependencies(expr: PyExpr) -> PyResult<PyPertMultichain> {
    let inner = expr.inner().clone();

    let e_ref = inner.as_any().downcast_ref::<TwoElecEnergy>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "two_elec_energy_dependencies() expected a TwoElecEnergy expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyPertMultichain::new(e_ref.dependencies().clone()))
}

/// Return derivative for a TwoElecEnergy.
///
/// Errors if the input expression is not a TwoElecEnergy.
#[pyfunction]
pub fn two_elec_energy_derivative(expr: PyExpr) -> PyResult<PyPertMultichain> {
    let inner = expr.inner().clone();

    let e_ref = inner.as_any().downcast_ref::<TwoElecEnergy>().ok_or_else(|| {
        to_pyerr(TinnedError::ExpressionError {
            message: "two_elec_energy_derivative() expected a TwoElecEnergy expression",
            expression: inner.to_string(),
            source: None,
        })
    })?;

    Ok(PyPertMultichain::new(e_ref.derivative().clone()))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(two_elec_energy_new, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_energy_from_operator, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_energy_name, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_energy_inner_density, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_energy_outer_density, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_energy_allow_density_swap, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_energy_dependencies, m)?)?;
    m.add_function(wrap_pyfunction!(two_elec_energy_derivative, m)?)?;

    Ok(())
}
