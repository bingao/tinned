use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::ZeroOperator;
use crate::public::{differentiate_expr, downcast_from_arc, downcast_from_ref, unreachable_error};

impl_nullary_expr_type!(LagMultiplier, LagMultiplierBuilder, false, false);
impl_nullary_expr_traits!(LagMultiplier, false, false);
