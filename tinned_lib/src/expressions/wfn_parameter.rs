use crate::core::Expr;

impl_nullary_expr_type!(WfnParameter, WfnParameterBuilder, false, false);
impl_nullary_expr_traits!(WfnParameter, false, false);

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::perturbations::pert_multichain::test_utils::make_pert_multichain;

    pub const TEST_OPER_NAME: &str = "psi";

    #[inline]
    pub fn make_wfn_parameter(name: impl Into<String>) -> std::sync::Arc<dyn Expr> {
        let name: String = name.into();
        if name.is_empty() {
            let deriv = make_pert_multichain(2u32, 10u32, 1u32, 10u32);
            WfnParameter::builder(random_alphanumeric(TEST_OPER_NAME.len() as u32 + 1))
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

    test_nullary_expr!(WfnParameter, TEST_OPER_NAME, make_wfn_parameter, false, false);
}
