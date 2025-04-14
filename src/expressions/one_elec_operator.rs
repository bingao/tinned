use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::perturbations::{PertMultichain, Perturbation};

impl_nullary_oper_type!(OneElecOperator, OneElecOperatorBuilder, true, false);
impl_nullary_oper_traits!(OneElecOperator, true, false);

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;
    use crate::expressions::{Number, Symbol};
    use crate::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_nullary_oper!(OneElecOperator, true, false);
}
