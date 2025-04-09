/// Represents exchange-correlation energy in a grid-based formulation.
/// xc_energy represents XC energy or its derivatives evaluated at grid points.
use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::expressions::{Add, Mul};
use crate::perturbations::{pert_multichain_hash_key, PertMultichain, Perturbation};
use crate::utils::{build_xc_density, downcast_expr, fmt_mul, intern, validate_xc_inputs};

impl_exch_corr_type!(ExchCorrEnergy, ExchCorrEnergyBuilder, xc_energy, true);
impl_exch_corr_traits!(ExchCorrEnergy, xc_energy, Mul, Add, fmt_mul, true);
