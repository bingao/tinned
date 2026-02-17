use pyo3::prelude::*;

pub mod number_tolerance;

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    number_tolerance::register(m)?;
    Ok(())
}
