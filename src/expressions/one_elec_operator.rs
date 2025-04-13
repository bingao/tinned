use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::perturbations::{PertMultichain, Perturbation};

impl_nullary_oper_type!(OneElecOperator, OneElecOperatorBuilder, true, false);
impl_nullary_oper_traits!(OneElecOperator, true, false);
