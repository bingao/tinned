use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::perturbations::{
    pert_multichain_display, pert_multichain_hash_key, PertMultichain, Perturbation,
};
use crate::utils::intern;

impl_nullary_oper_type!(LagMultiplier, LagMultiplierBuilder, false, false);
impl_nullary_oper_traits!(LagMultiplier, false, false);
