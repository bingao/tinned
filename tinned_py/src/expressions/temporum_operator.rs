use pyo3::prelude::*;

use tinned::TemporumOperator;

use crate::core::{errors::to_pyerr, expr::PyExpr};
use crate::perturbations::pert_multichain::PyPertMultichain;

/// Create a TemporumOperator expression.
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
pub fn temporum_operator_new(argument: PyExpr, is_forward: Option<bool>) -> PyResult<PyExpr> {
    let mut b = TemporumOperator::builder(argument.inner().clone());

    if let Some(v) = is_forward {
        b = b.is_forward(v);
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = temporum_operator_is_forward,
    fn_doc = impl_expr_getter_doc!("i*d/dt (True), or -i*d/dt (False)", TemporumOperator),
    expr_ty = TemporumOperator,
    out_ty = bool,
    body = |op: &TemporumOperator| Ok(op.is_forward())
);

impl_expr_getter_interface!(
    fn_name = temporum_operator_argument,
    fn_doc = impl_expr_getter_doc!("argument", TemporumOperator),
    expr_ty = TemporumOperator,
    out_ty = PyExpr,
    body = |op: &TemporumOperator| Ok(PyExpr::new(op.argument().clone()))
);

impl_expr_getter_interface!(
    fn_name = temporum_operator_derivative,
    fn_doc = impl_expr_getter_doc!("derivative", TemporumOperator),
    expr_ty = TemporumOperator,
    out_ty = PyPertMultichain,
    body = |op: &TemporumOperator| {
        let derivative = op.derivative().map_err(to_pyerr)?;
        Ok(PyPertMultichain::new(derivative.clone()))
    }
);

impl_expr_getter_interface!(
    fn_name = temporum_operator_frequency,
    fn_doc = impl_expr_getter_doc!("frequency", TemporumOperator),
    expr_ty = TemporumOperator,
    out_ty = PyExpr,
    body = |op: &TemporumOperator| {
        let freq = op.frequency().map_err(to_pyerr)?;
        Ok(PyExpr::new(freq))
    }
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(temporum_operator_new, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_operator_is_forward, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_operator_argument, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_operator_derivative, m)?)?;
    m.add_function(wrap_pyfunction!(temporum_operator_frequency, m)?)?;

    Ok(())
}
