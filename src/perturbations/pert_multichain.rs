use std::sync::Arc;

use crate::perturbations::Perturbation;

/// Perturbation multichain: unique perturbations with associated differentiation orders.
pub type PertMultichain = std::collections::BTreeMap<Arc<Perturbation>, u32>;

/// Generates a compact string suitable for hashing a perturbation multichain.
#[inline]
pub fn pert_multichain_hash_key(multichain: &PertMultichain) -> String {
    multichain
        .iter()
        .map(|(p, order)| format!("{}^{}", p.hash_key(), order))
        .collect::<Vec<_>>()
        .join(",")
}

/// Produces a human-readable string for a perturbation multichain.
#[inline]
pub fn pert_multichain_display(multichain: &PertMultichain) -> String {
    multichain.iter().map(|(p, order)| format!("{}^{}", p, order)).collect::<Vec<_>>().join(", ")
}

/// Returns the order of a perturbation in the multichain, or 0 if absent.
#[inline]
pub fn get_perturbation_order(multichain: &PertMultichain, p: &Arc<Perturbation>) -> u32 {
    *multichain.get(p).unwrap_or(&0)
}

/// Returns true if `chain` is a sub-chain of `super_chain`.
/// That is, all perturbations in `chain` appear in `super_chain`
/// with at least the same order.
#[inline]
pub fn is_sub_multichain(chain: &PertMultichain, super_chain: &PertMultichain) -> bool {
    chain.iter().all(|(p, &order)| order <= get_perturbation_order(super_chain, p))
}
