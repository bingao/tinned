use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Number, ZeroOperator};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{downcast_from_arc, downcast_from_ref};

impl_nullary_expr_type!(WfnParameter, WfnParameterBuilder, false, false);
impl_nullary_expr_traits!(WfnParameter, false, false);

#[cfg(test)]
const DEFAULT_OPER_NAME: &str = "psi";

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::perturbations::pert_multichain::test_utils::make_pert_multichain;

    #[inline]
    pub fn make_wfn_parameter(name: impl Into<String>) -> Arc<dyn Expr> {
        let name: String = name.into();
        if name.is_empty() {
            let deriv = make_pert_multichain(2u32, 10u32, 1u32, 10u32);
            WfnParameter::builder(random_alphanumeric(DEFAULT_OPER_NAME.len() as u32 + 1))
                .derivative(deriv)
                .build()
                .unwrap()
        } else {
            let deriv = make_pert_multichain(0u32, 0u32, 1u32, 0u32);
            WfnParameter::builder(name).derivative(deriv).build().unwrap()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::perturbations::pert_multichain::test_utils::make_pert_multichain;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_nullary_expr!(WfnParameter, DEFAULT_OPER_NAME, make_wfn_parameter, false, false);
}
