use pyo3::prelude::*;

use tinned::Transpose;

use crate::core::{errors::to_pyerr, expr::PyExpr};

/// Create a Transpose expression from an argument.
///
/// Args:
///   argument: Argument expression (non-scalar).
///   is_hermitian: Optional Hermitian transpose flag (default False).
///
/// Returns:
///   A PyExpr wrapping the constructed Transpose (interned).
///
/// Notes:
///   - argument must be non-scalar.
///   - If argument is ZeroOperator, the result is ZeroOperator.
#[pyfunction]
#[pyo3(signature = (
    argument,
    is_hermitian=false
))]
pub fn transpose_new(argument: &PyExpr, is_hermitian: bool) -> PyResult<PyExpr> {
    let rust_argument = argument.inner().clone();
    let out = Transpose::new(rust_argument, is_hermitian).map_err(to_pyerr)?;

    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = transpose_argument,
    fn_doc = impl_expr_getter_doc!("argument", Transpose),
    expr_ty = Transpose,
    out_ty = PyExpr,
    body = |op: &Transpose| Ok(PyExpr::new(op.argument().clone()))
);

impl_expr_getter_interface!(
    fn_name = transpose_is_hermitian,
    fn_doc = impl_expr_getter_doc!(
        "standard transpose (False) or Hermitian transpose (True)",
        Transpose
    ),
    expr_ty = Transpose,
    out_ty = bool,
    body = |op: &Transpose| Ok(op.is_hermitian())
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(transpose_new, m)?)?;
    m.add_function(wrap_pyfunction!(transpose_argument, m)?)?;
    m.add_function(wrap_pyfunction!(transpose_is_hermitian, m)?)?;
    Ok(())
}
