use paste::paste;
use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::ExchCorrEnergy;
use tinned::public::generic_error;

use crate::c_support::{
    tinned_string_from_cstr, with_downcast_cstr, with_downcast_expr, with_downcast_pert_multichain,
};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::PertMultichainBox;

impl_exch_corr_ffi!(ExchCorrEnergy, exch_corr_energy, xc_energy);
