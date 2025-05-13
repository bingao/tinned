use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::ZeroOperator;
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{differentiate_expr, downcast_from_arc, downcast_from_ref};

impl_nullary_expr_type!(LagMultiplier, LagMultiplierBuilder, false, false);
impl_nullary_expr_traits!(LagMultiplier, false, false);
