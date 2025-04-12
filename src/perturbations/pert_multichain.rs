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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::{Number, Symbol};

    #[test]
    fn test_pert_multichain_behavior() {
        let p1 = Perturbation::new("p1", Number::Integer(1).into());
        let p2 = Perturbation::new("p2", Symbol::new("x"));

        let mut chain = PertMultichain::new();
        chain.insert(p1.clone(), 1);
        chain.insert(p2.clone(), 2);

        // Order is stored correctly
        assert_eq!(get_perturbation_order(&chain, &p1), 1);
        assert_eq!(get_perturbation_order(&chain, &p2), 2);

        // Absent perturbation returns 0
        let absent = Perturbation::new("absent", Number::Real(3.14).into());
        assert_eq!(get_perturbation_order(&chain, &absent), 0);

        // Display and hash key format
        let display = pert_multichain_display(&chain);
        assert!(display.contains("p1^1"));
        assert!(display.contains("p2^2"));

        let hash_key = pert_multichain_hash_key(&chain);
        assert!(hash_key.contains("p1^1"));
        assert!(hash_key.contains("p2^2"));
    }

    #[test]
    fn test_is_sub_multichain() {
        let p1 = Perturbation::new("p1", Number::Integer(1).into());
        let p2 = Perturbation::new("p2", Symbol::new("x"));

        let mut super_chain = PertMultichain::new();
        super_chain.insert(p1.clone(), 2);
        super_chain.insert(p2.clone(), 3);

        let mut sub_chain = PertMultichain::new();
        sub_chain.insert(p1.clone(), 1);

        let mut exact_chain = PertMultichain::new();
        exact_chain.insert(p1.clone(), 2);
        exact_chain.insert(p2.clone(), 3);

        let mut invalid_sub_chain = PertMultichain::new();
        invalid_sub_chain.insert(p2.clone(), 4); // Too high order

        assert!(is_sub_multichain(&sub_chain, &super_chain));
        assert!(is_sub_multichain(&exact_chain, &super_chain));
        assert!(!is_sub_multichain(&invalid_sub_chain, &super_chain));
    }
}
