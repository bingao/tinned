use pyo3::prelude::*;

use tinned::Symbol;

use crate::core::expr::PyExpr;

/// Build a Symbol from a name.
///
/// Args:
///   name: Symbol name.
///
/// Returns:
///   A PyExpr wrapping the constructed symbol (interned).
#[pyfunction]
pub fn symbol_new(name: String) -> PyResult<PyExpr> {
    let out = Symbol::new(name);
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = symbol_name,
    fn_doc = impl_expr_getter_doc!("name", Symbol),
    expr_ty = Symbol,
    out_ty = String,
    body = |op: &Symbol| Ok(op.name().to_string())
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(symbol_new, m)?)?;
    m.add_function(wrap_pyfunction!(symbol_name, m)?)?;

    Ok(())
}
