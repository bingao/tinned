use pyo3::prelude::*;

pub mod add;
pub mod adjoint_map;
pub mod composition;
pub mod conjugate;
pub mod dot_product;
pub mod exch_corr_energy;
pub mod exch_corr_potential;
pub mod exp_adjoint_map;
pub mod hermitian_transpose;
pub mod lag_multiplier;
pub mod matrix_add;
pub mod matrix_mul;
pub mod mul;
pub mod non_elec_function;
pub mod number;
pub mod one_elec_operator;
pub mod power;
pub mod residue_parameter;
pub mod symbol;
pub mod temporum_operator;
pub mod temporum_overlap;
pub mod trace;
pub mod transpose;
pub mod two_elec_energy;
pub mod two_elec_operator;
pub mod wfn_parameter;
pub mod zero_operator;

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    add::register(m)?;
    adjoint_map::register(m)?;
    composition::register(m)?;
    conjugate::register(m)?;
    dot_product::register(m)?;
    exch_corr_energy::register(m)?;
    exch_corr_potential::register(m)?;
    exp_adjoint_map::register(m)?;
    hermitian_transpose::register(m)?;
    lag_multiplier::register(m)?;
    matrix_add::register(m)?;
    matrix_mul::register(m)?;
    mul::register(m)?;
    non_elec_function::register(m)?;
    number::register(m)?;
    one_elec_operator::register(m)?;
    power::register(m)?;
    residue_parameter::register(m)?;
    symbol::register(m)?;
    temporum_operator::register(m)?;
    temporum_overlap::register(m)?;
    trace::register(m)?;
    transpose::register(m)?;
    two_elec_energy::register(m)?;
    two_elec_operator::register(m)?;
    wfn_parameter::register(m)?;
    zero_operator::register(m)?;

    Ok(())
}
