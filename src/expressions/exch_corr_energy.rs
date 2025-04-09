/// Represents exchange-correlation energy in a grid-based formulation.
/// xc_energy represents XC energy or its derivatives evaluated at grid points.
use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::perturbations::{pert_multichain_hash_key, PertMultichain, Perturbation};
use crate::utils::{
    build_xc_density, downcast_from_arc, downcast_from_ref, fmt_mul, intern, validate_xc_inputs,
};

impl_exch_corr_type!(ExchCorrEnergy, ExchCorrEnergyBuilder, xc_energy, true);
impl_exch_corr_traits!(
    ExchCorrEnergy,
    xc_energy,
    crate::expressions::Mul,
    crate::expressions::Add,
    fmt_mul,
    true
);
