use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::perturbations::{PertMultichain, Perturbation};

impl_nullary_oper_type!(LagMultiplier, LagMultiplierBuilder, false, false);
impl_nullary_oper_traits!(LagMultiplier, false, false);
