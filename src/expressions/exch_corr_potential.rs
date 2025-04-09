/// Represents exchange-correlation potential in a grid-based formulation.
/// xc_potential represents XC potential or its derivatives evaluated at grid points.
use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::perturbations::{pert_multichain_hash_key, PertMultichain, Perturbation};
use crate::utils::{
    build_xc_density, downcast_from_arc, downcast_from_ref, fmt_matrix_mul, intern,
    validate_xc_inputs,
};

impl_exch_corr_type!(ExchCorrPotential, ExchCorrPotentialBuilder, xc_potential, false);
impl_exch_corr_traits!(
    ExchCorrPotential,
    xc_potential,
    crate::expressions::MatrixMul,
    crate::expressions::MatrixAdd,
    fmt_matrix_mul,
    false
);
