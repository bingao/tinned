use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Number, ZeroOperator};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{downcast_from_arc, downcast_from_ref};

impl_nullary_expr_type!(LagMultiplier, LagMultiplierBuilder, false, false);
impl_nullary_expr_traits!(LagMultiplier, false, false);
