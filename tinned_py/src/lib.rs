use pyo3::prelude::*;

#[macro_use]
mod macros;

mod core;
mod expressions;
mod perturbations;
mod public;

#[pymodule]
fn tinned(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Register exceptions and base Expr
    core::register(m)?;

    // Register perturbations and perturbation multichains
    perturbations::register(m)?;

    // Create and attach tinned.expressions
    let expressions_mod = PyModule::new(m.py(), "expressions")?;
    expressions::register(&expressions_mod)?;
    m.add_submodule(&expressions_mod)?;

    // Public classes and functions
    public::register(m)?;

    Ok(())
}
