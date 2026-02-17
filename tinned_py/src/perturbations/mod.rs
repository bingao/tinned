use pyo3::prelude::*;

pub mod pert_multichain;
pub mod perturbation;

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    perturbation::register(m)?;
    pert_multichain::register(m)?;
    Ok(())
}
