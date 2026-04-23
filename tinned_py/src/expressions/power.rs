use pyo3::prelude::*;

use tinned::Power;

use crate::core::{errors::to_pyerr, expr::PyExpr};

/// Build a power expression from base and exponent.
///
/// Args:
///   base: Scalar base expression.
///   exponent: Integer exponent.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
pub fn power_new(base: &PyExpr, exponent: i64) -> PyResult<PyExpr> {
    let base_inner = base.inner().clone();

    let out = Power::new(base_inner, exponent).map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = power_base,
    fn_doc = impl_expr_getter_doc!("base", Power),
    expr_ty = Power,
    out_ty = PyExpr,
    body = |power: &Power| Ok(PyExpr::new(power.base().clone()))
);

impl_expr_getter_interface!(
    fn_name = power_exponent,
    fn_doc = impl_expr_getter_doc!("exponent", Power),
    expr_ty = Power,
    out_ty = i64,
    body = |power: &Power| Ok(power.exponent())
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(power_new, m)?)?;
    m.add_function(wrap_pyfunction!(power_base, m)?)?;
    m.add_function(wrap_pyfunction!(power_exponent, m)?)?;

    Ok(())
}
