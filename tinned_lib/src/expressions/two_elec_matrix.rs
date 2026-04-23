use crate::core::Expr;

// Two-electron second-quantization matrix, see (1.4.37) of the pink bible
// "Molecular Electronic-Structure Theory", which can be contracted with
// two-electron density matrix (1.7.6)
impl_nullary_expr_type!(TwoElecMatrix, TwoElecMatrixBuilder, true, false);
impl_nullary_expr_traits!(TwoElecMatrix, true, false);

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::perturbations::pert_multichain::test_utils::{
        make_pert_multichain, make_super_multichain,
    };

    pub const TEST_OPER_NAME: &str = "g";

    #[inline]
    pub fn make_two_elec_matrix(name: impl Into<String>) -> std::sync::Arc<dyn Expr> {
        let name: String = name.into();
        if name.is_empty() {
            let deriv = make_pert_multichain(2u32, 10u32, 1u32, 10u32);
            let deps = make_super_multichain(&deriv, 1u32);
            TwoElecMatrix::builder(random_alphanumeric(TEST_OPER_NAME.len() as u32 + 1))
                .is_perturbing(false)
                .dependencies(deps)
                .derivative(deriv)
                .build()
                .unwrap()
        } else {
            let deriv = make_pert_multichain(0u32, 0u32, 1u32, 0u32);
            let deps = make_super_multichain(&deriv, 1u32);
            TwoElecMatrix::builder(name)
                .is_perturbing(false)
                .dependencies(deps)
                .derivative(deriv)
                .build()
                .unwrap()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;

    test_nullary_expr!(TwoElecMatrix, TEST_OPER_NAME, make_two_elec_matrix, true, false);
}
