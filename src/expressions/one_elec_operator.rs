use std::collections::HashSet;
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::ZeroOperator;
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{differentiate_expr, downcast_from_arc, downcast_from_ref, unreachable_error};

impl_nullary_expr_type!(OneElecOperator, OneElecOperatorBuilder, true, false);
impl_nullary_expr_traits!(OneElecOperator, true, false);

#[cfg(test)]
const DEFAULT_OPER_NAME: &str = "op(1el)";

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::perturbations::pert_multichain::test_utils::{
        make_pert_multichain, make_super_multichain,
    };

    #[inline]
    pub fn make_one_elec_operator(name: impl Into<String>) -> Arc<dyn Expr> {
        let name: String = name.into();
        if name.is_empty() {
            let deriv = make_pert_multichain(2u32, 10u32, 1u32, 10u32);
            let deps = make_super_multichain(&deriv, 1u32);
            OneElecOperator::builder(random_alphanumeric(DEFAULT_OPER_NAME.len() as u32 + 1))
                .dependencies(deps)
                .derivative(deriv)
                .build()
                .unwrap()
        } else {
            let deriv = make_pert_multichain(0u32, 0u32, 1u32, 0u32);
            let deps = make_super_multichain(&deriv, 1u32);
            OneElecOperator::builder(name).dependencies(deps).derivative(deriv).build().unwrap()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::perturbations::pert_multichain::test_utils::{
        make_pert_multichain, make_super_multichain,
    };
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_nullary_expr!(OneElecOperator, DEFAULT_OPER_NAME, make_one_elec_operator, true, false);
}
