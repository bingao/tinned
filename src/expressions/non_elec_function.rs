use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::perturbations::{PertMultichain, Perturbation};

impl_nullary_oper_type!(NonElecFunction, NonElecFunctionBuilder, true, true);
impl_nullary_oper_traits!(NonElecFunction, true, true);

#[cfg(test)]
const DEFAULT_FUN_NAME: &str = "nel";

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::perturbations::pert_multichain::test_utils::{
        make_pert_multichain, make_super_multichain,
    };

    #[inline]
    pub fn make_non_elec_function(name: impl Into<String>) -> Arc<dyn Expr> {
        let name: String = name.into();
        if name.is_empty() {
            let deriv = make_pert_multichain(2u32, 10u32, 1u32, 10u32);
            let deps = make_super_multichain(&deriv, 1u32);
            NonElecFunction::builder(random_alphanumeric(DEFAULT_FUN_NAME.len() as u32 + 1))
                .dependencies(deps)
                .derivative(deriv)
                .build()
                .unwrap()
        } else {
            let deriv = make_pert_multichain(0u32, 0u32, 1u32, 0u32);
            let deps = make_super_multichain(&deriv, 1u32);
            NonElecFunction::builder(name).dependencies(deps).derivative(deriv).build().unwrap()
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
    use crate::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_nullary_oper!(NonElecFunction, DEFAULT_FUN_NAME, make_non_elec_function, true, true);
}
