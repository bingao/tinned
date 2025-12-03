use paste::paste;
use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::ExchCorrEnergy;
use tinned::public::generic_error;

use crate::c_support::{ffi_map_expr_as, tinned_string_from_cstr, tinned_string_to_cstr};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::{PertMultichainBox, PertMultichainHandle};

impl_exch_corr_ffi!(ExchCorrEnergy, exch_corr_energy, xc_energy);
