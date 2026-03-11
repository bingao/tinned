use crate::core::Expr;
use crate::core::expr_internal::sealed::ExprInternal;

impl_nullary_expr_type!(OneElecOperator, OneElecOperatorBuilder, true, false);
impl_nullary_expr_traits!(OneElecOperator, true, false);

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::perturbations::pert_multichain::test_utils::{
        make_pert_multichain, make_super_multichain,
    };

    pub const TEST_OPER_NAME: &str = "op(1el)";

    #[inline]
    pub fn make_one_elec_operator(name: impl Into<String>) -> std::sync::Arc<dyn Expr> {
        let name: String = name.into();
        if name.is_empty() {
            let deriv = make_pert_multichain(2u32, 10u32, 1u32, 10u32);
            let deps = make_super_multichain(&deriv, 1u32);
            OneElecOperator::builder(random_alphanumeric(TEST_OPER_NAME.len() as u32 + 1))
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

    test_nullary_expr!(OneElecOperator, TEST_OPER_NAME, make_one_elec_operator, true, false);
}
