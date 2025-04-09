use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::expressions::ZeroOperator;
use crate::perturbations::{
    is_sub_multichain, pert_multichain_display, pert_multichain_hash_key, PertMultichain,
    Perturbation,
};
use crate::utils::intern;

impl_nullary_oper_type!(OneElecOperator, OneElecOperatorBuilder, true, false);
impl_nullary_oper_traits!(OneElecOperator, true, false);
