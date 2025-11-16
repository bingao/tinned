use paste::paste;
use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::ExchCorrPotential;
use tinned::public::generic_error;

use crate::c_support::{
    tinned_string_from_cstr, ffi_map_expr_as, tinned_string_to_cstr,
};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::{PertMultichainHandle, PertMultichainBox};

impl_exch_corr_ffi!(ExchCorrPotential, exch_corr_potential, xc_potential);
