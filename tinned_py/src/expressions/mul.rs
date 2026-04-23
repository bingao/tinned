use pyo3::prelude::*;
use std::vec::Vec;

use tinned::Mul;

use crate::core::{errors::to_pyerr, expr::PyExpr};

/// Build a multiplication from scalar terms.
///
/// Args:
///   terms: Scalar terms.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
pub fn mul_new(terms: Vec<PyExpr>) -> PyResult<PyExpr> {
    let rust_terms = terms.into_iter().map(|t| t.inner().clone()).collect();

    let out = Mul::new(rust_terms).map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = mul_coefficient,
    fn_doc = impl_expr_getter_doc!("coefficient", Mul),
    expr_ty = Mul,
    out_ty = PyExpr,
    body = |mul: &Mul| {
        // Convert &Number to Arc<dyn Expr> and wrap in PyExpr.
        let coeff_expr = mul.coefficient().clone().into();

        Ok(PyExpr::new(coeff_expr))
    }
);

impl_expr_getter_interface!(
    fn_name = mul_factors,
    fn_doc = impl_expr_getter_doc!("factors (excluding the coefficient)", Mul),
    expr_ty = Mul,
    out_ty = Vec<PyExpr>,
    body = |mul: &Mul| Ok(mul.factors().iter().cloned().map(PyExpr::new).collect())
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(mul_new, m)?)?;
    m.add_function(wrap_pyfunction!(mul_coefficient, m)?)?;
    m.add_function(wrap_pyfunction!(mul_factors, m)?)?;

    Ok(())
}
