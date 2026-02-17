use pyo3::prelude::*;

pub mod errors;
pub mod expr;

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    errors::register(m)?;
    expr::register(m)?;
    Ok(())
}
