use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::perturbations::{PertMultichain, Perturbation};

impl_nullary_oper_type!(NonElecFunction, NonElecFunctionBuilder, true, true);
impl_nullary_oper_traits!(NonElecFunction, true, true);
