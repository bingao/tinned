use crate::core::Expr;
use crate::core::expr_internal::sealed::ExprInternal;

impl_nullary_expr_type!(LagMultiplier, LagMultiplierBuilder, false, false);
impl_nullary_expr_traits!(LagMultiplier, false, false);
