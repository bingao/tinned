use pyo3::prelude::*;

use tinned::TimeEvolution;

use crate::core::{errors::to_pyerr, expr::PyExpr};
use crate::perturbations::pert_multichain::PyPertMultichain;

/// Create a TimeEvolution expression.
///
/// Args:
///   argument: Non-scalar argument expression.
///   is_forward: Optional bool. If True uses i*d/dt, otherwise -i*d/dt.
///               Defaults to True.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
#[pyo3(signature = (argument, is_forward=None))]
pub fn time_evolution_new(argument: PyExpr, is_forward: Option<bool>) -> PyResult<PyExpr> {
    let mut b = TimeEvolution::builder(argument.inner().clone());

    if let Some(v) = is_forward {
        b = b.is_forward(v);
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = time_evolution_is_forward,
    fn_doc = impl_expr_getter_doc!("i*d/dt (True), or -i*d/dt (False)", TimeEvolution),
    expr_ty = TimeEvolution,
    out_ty = bool,
    body = |op: &TimeEvolution| Ok(op.is_forward())
);

impl_expr_getter_interface!(
    fn_name = time_evolution_argument,
    fn_doc = impl_expr_getter_doc!("argument", TimeEvolution),
    expr_ty = TimeEvolution,
    out_ty = PyExpr,
    body = |op: &TimeEvolution| Ok(PyExpr::new(op.argument().clone()))
);

impl_expr_getter_interface!(
    fn_name = time_evolution_derivative,
    fn_doc = impl_expr_getter_doc!("derivative", TimeEvolution),
    expr_ty = TimeEvolution,
    out_ty = PyPertMultichain,
    body = |op: &TimeEvolution| {
        let derivative = op.derivative().map_err(to_pyerr)?;
        Ok(PyPertMultichain::new(derivative.clone()))
    }
);

impl_expr_getter_interface!(
    fn_name = time_evolution_frequency,
    fn_doc = impl_expr_getter_doc!("frequency", TimeEvolution),
    expr_ty = TimeEvolution,
    out_ty = PyExpr,
    body = |op: &TimeEvolution| {
        let freq = op.frequency().map_err(to_pyerr)?;
        Ok(PyExpr::new(freq))
    }
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(time_evolution_new, m)?)?;
    m.add_function(wrap_pyfunction!(time_evolution_is_forward, m)?)?;
    m.add_function(wrap_pyfunction!(time_evolution_argument, m)?)?;
    m.add_function(wrap_pyfunction!(time_evolution_derivative, m)?)?;
    m.add_function(wrap_pyfunction!(time_evolution_frequency, m)?)?;

    Ok(())
}
