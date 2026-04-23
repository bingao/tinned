use crate::core::Expr;

impl_nullary_expr_type!(NonElecFunction, NonElecFunctionBuilder, true, true);
impl_nullary_expr_traits!(NonElecFunction, true, true);

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::perturbations::pert_multichain::test_utils::{
        make_pert_multichain, make_super_multichain,
    };

    pub const TEST_FUN_NAME: &str = "nel";

    #[inline]
    pub fn make_non_elec_function(name: impl Into<String>) -> std::sync::Arc<dyn Expr> {
        let name: String = name.into();
        if name.is_empty() {
            let deriv = make_pert_multichain(2u32, 10u32, 1u32, 10u32);
            let deps = make_super_multichain(&deriv, 1u32);
            NonElecFunction::builder(random_alphanumeric(TEST_FUN_NAME.len() as u32 + 1))
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

    test_nullary_expr!(NonElecFunction, TEST_FUN_NAME, make_non_elec_function, true, true);
}
