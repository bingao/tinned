use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::perturbations::Perturbation;

/// Perturbation multichain: unique perturbations with associated differentiation orders.
#[derive(Clone, Debug)]
pub struct PertMultichain(Arc<Mutex<BTreeMap<Arc<Perturbation>, u32>>>);

impl PertMultichain {
    /// Creates a new and empty perturbation multichain.
    #[inline]
    pub fn new() -> Self {
        Self(Arc::new(Mutex::new(BTreeMap::new())))
    }

    /// Creates a new perturbation multichain from a given `BTreeMap`.
    #[inline]
    pub fn from_map(map: BTreeMap<Arc<Perturbation>, u32>) -> Self {
        PertMultichain(Arc::new(Mutex::new(map)))
    }

    /// Returns true if the perturbation multichain is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.lock().unwrap().is_empty()
    }

    /// Returns the order of a perturbation in the multichain, or 0 if absent.
    #[inline]
    pub fn get_order(&self, p: &Arc<Perturbation>) -> u32 {
        let map = self.0.lock().unwrap();
        *map.get(p).unwrap_or(&0)
    }

    /// Returns all perturbations in the multichain and meanwhile preserves the order.
    pub fn keys(&self) -> Vec<Arc<Perturbation>> {
        let map = self.0.lock().unwrap();
        map.keys().cloned().collect()
    }

    /// Returns a cloned internal map of the perturbation multichain.
    #[inline]
    pub fn get_map_clone(&self) -> BTreeMap<Arc<Perturbation>, u32> {
        self.0.lock().unwrap().clone()
    }

    /// Inserts a perturbation into the multichain, or increases the order by 1
    /// if the perturbation already exists in the multichain.
    #[inline]
    pub fn insert(&mut self, p: &Arc<Perturbation>) {
        let mut map = self.0.lock().unwrap();
        let entry = map.entry(p.clone()).or_insert(0);
        *entry += 1;
    }

    /// Returns true if `chain` is a super-multichain of `subchain`. That is,
    /// all perturbations in `subchain` appear in `chain` with at least the
    /// same order.
    #[inline]
    pub fn is_subchain(&self, subchain: &PertMultichain) -> bool {
        let map = self.0.lock().unwrap().clone();
        let submap = subchain.0.lock().unwrap().clone();

        submap.iter().all(|(p, &order)| map.get(p).copied().unwrap_or(0) >= order)
    }

    /// Returns true if `chain` is a sub-multichain of `superchain`. That is,
    /// all perturbations in `chain` appear in `superchain` with at least the
    /// same order.
    #[inline]
    pub fn is_superchain(&self, superchain: &PertMultichain) -> bool {
        let map = self.0.lock().unwrap().clone();
        let supermap = superchain.0.lock().unwrap().clone();

        map.iter().all(|(p, &order)| supermap.get(p).copied().unwrap_or(0) >= order)
    }

    /// Generates a compact string suitable for hashing a perturbation multichain.
    #[inline]
    pub fn hash_key(&self) -> String {
        let map = self.0.lock().unwrap();
        let mut parts = Vec::new();
        for (pert, order) in map.iter() {
            parts.push(format!("{}^{}", pert.name(), order));
        }
        parts.join(",")
    }
}

impl PartialEq for PertMultichain {
    fn eq(&self, other: &Self) -> bool {
        // Clone to avoid double-locks, nested locks and deadlock
        let self_map = self.0.lock().unwrap().clone();
        let other_map = other.0.lock().unwrap().clone();
        self_map == other_map
    }
}

impl Eq for PertMultichain {}

impl std::fmt::Display for PertMultichain {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let map = self.0.lock().unwrap();
        let parts: Vec<String> =
            map.iter().map(|(pert, order)| format!("{}^{}", pert.name(), order)).collect();
        f.write_str(&parts.join(","))
    }
}

// Intermediate format for serialization
#[derive(Serialize, Deserialize)]
struct PertEntry {
    perturbation: Arc<Perturbation>,
    order: u32,
}

impl Serialize for PertMultichain {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let map = self.0.lock().unwrap();
        let entries: Vec<PertEntry> = map
            .iter()
            .map(|(p, &order)| PertEntry {
                perturbation: Arc::clone(p),
                order,
            })
            .collect();
        entries.serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for PertMultichain {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let entries: Vec<PertEntry> = Vec::deserialize(deserializer)?;
        let mut map = BTreeMap::new();

        for entry in entries {
            map.insert(entry.perturbation, entry.order);
        }

        Ok(PertMultichain(Arc::new(Mutex::new(map))))
    }
}

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::perturbations::perturbation::test_utils::*;

    #[inline]
    pub fn make_pert_multichain(
        len_name: u32,
        val_range: u32,
        min_order: u32,
        order_range: u32,
    ) -> PertMultichain {
        let mut map: BTreeMap<Arc<Perturbation>, u32> = std::collections::BTreeMap::new();

        let pert: Vec<Arc<Perturbation>> = vec![
            make_perturbation_i64(len_name, val_range),
            make_perturbation_f64(len_name, val_range),
            make_perturbation_complex(len_name, val_range),
            make_perturbation_rational(len_name, val_range),
            make_perturbation_symbol(len_name, val_range),
        ];

        for p in &pert {
            let order: u32 = rand::random_range(min_order..=min_order + order_range);
            map.insert(Arc::clone(p), order);
        }

        PertMultichain::from_map(map)
    }

    #[inline]
    pub fn make_super_multichain(subchain: &PertMultichain, min_increase: u32) -> PertMultichain {
        let increment = if min_increase == 0 {
            1
        } else {
            min_increase
        };

        let mut supermap = BTreeMap::new();
        let map = subchain.get_map_clone();

        for (pert, order) in map.iter() {
            supermap.insert(Arc::clone(pert), order + increment);
        }

        PertMultichain::from_map(supermap)
    }
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;
    use crate::perturbations::perturbation::test_utils::*;

    test_struct_safety!(PertMultichain);

    #[test]
    fn test_setter_and_getter() {
        let mut chain = PertMultichain::new();

        let p1 = make_perturbation_i64(2u32, 10u32);

        assert_eq!(chain.get_order(&p1), 0);

        chain.insert(&p1);
        assert_eq!(chain.get_order(&p1), 1);

        chain.insert(&p1);
        assert_eq!(chain.get_order(&p1), 2);

        let mut keys: Vec<Arc<Perturbation>> = chain.keys();

        assert_eq!(keys, vec![p1.clone()]);

        let p2 = make_perturbation_f64(2u32, 10u32);
        let p3 = make_perturbation_complex(2u32, 10u32);
        let p4 = make_perturbation_rational(2u32, 10u32);
        let p5 = make_perturbation_symbol(2u32, 4u32);

        chain.insert(&p2);
        chain.insert(&p3);
        chain.insert(&p4);
        chain.insert(&p5);

        keys = chain.keys();

        let mut expected_keys: Vec<Arc<Perturbation>> = vec![p1, p2, p3, p4, p5];
        expected_keys.sort();

        assert_eq!(keys, expected_keys);
    }

    #[test]
    fn test_from_map() {
        let mut map = BTreeMap::new();

        let p1 = make_perturbation_i64(2u32, 10u32);
        let p2 = make_perturbation_f64(2u32, 10u32);
        let p3 = make_perturbation_complex(2u32, 10u32);
        let p4 = make_perturbation_rational(2u32, 10u32);
        let p5 = make_perturbation_symbol(2u32, 4u32);

        map.insert(p1.clone(), 2);
        map.insert(p2.clone(), 1);
        map.insert(p3.clone(), 3);
        map.insert(p4.clone(), 5);
        map.insert(p5.clone(), 4);

        let chain = PertMultichain::from_map(map);

        assert_eq!(chain.get_order(&p1), 2);
        assert_eq!(chain.get_order(&p2), 1);
        assert_eq!(chain.get_order(&p3), 3);
        assert_eq!(chain.get_order(&p4), 5);
        assert_eq!(chain.get_order(&p5), 4);
    }

    #[test]
    fn test_get_map_clone() {
        let mut chain = PertMultichain::new();

        let p1 = make_perturbation_i64(2u32, 10u32);
        let p2 = make_perturbation_f64(2u32, 10u32);
        let p3 = make_perturbation_complex(2u32, 10u32);
        let p4 = make_perturbation_rational(2u32, 10u32);
        let p5 = make_perturbation_symbol(2u32, 4u32);

        chain.insert(&p1);
        chain.insert(&p1);
        chain.insert(&p1);
        chain.insert(&p2);
        chain.insert(&p2);
        chain.insert(&p3);
        chain.insert(&p4);

        let map = chain.get_map_clone();

        assert_eq!(map.get(&p1), Some(&3));
        assert_eq!(map.get(&p2), Some(&2));
        assert_eq!(map.get(&p3), Some(&1));
        assert_eq!(map.get(&p4), Some(&1));
        assert_eq!(map.get(&p5), None);
    }

    #[test]
    fn test_chain_relationships() {
        let mut c1 = make_pert_multichain(2u32, 8u32, 0u32, 10u32);
        let mut c2 = make_super_multichain(&c1, 1u32);

        assert!(c1.is_superchain(&c2));
        assert!(!c1.is_subchain(&c2));
        assert!(!c2.is_superchain(&c1));
        assert!(c2.is_subchain(&c1));

        for p in c2.keys() {
            c1.insert(&p);
        }

        assert_eq!(c1, c2);

        assert!(c1.is_superchain(&c2));
        assert!(c1.is_subchain(&c2));
        assert!(c2.is_superchain(&c1));
        assert!(c2.is_subchain(&c1));

        c1.insert(&make_perturbation_i64(2u32, 10u32));
        c2.insert(&make_perturbation_f64(2u32, 10u32));

        assert!(!c1.is_superchain(&c2));
        assert!(!c1.is_subchain(&c2));
        assert!(!c2.is_superchain(&c1));
        assert!(!c2.is_subchain(&c1));
    }

    #[test]
    fn test_hash_key_and_display() {
        let mut map = BTreeMap::new();

        let p1 = make_perturbation_i64(2u32, 10u32);
        let p2 = make_perturbation_f64(2u32, 10u32);
        let p3 = make_perturbation_complex(2u32, 10u32);
        let p4 = make_perturbation_rational(2u32, 10u32);
        let p5 = make_perturbation_symbol(2u32, 4u32);

        map.insert(p1.clone(), 2);
        map.insert(p2.clone(), 1);
        map.insert(p3.clone(), 3);
        map.insert(p4.clone(), 5);
        map.insert(p5.clone(), 4);

        let chain = PertMultichain::from_map(map);

        let key = chain.hash_key();
        assert!(
            key.contains(&format!("{}^2", p1.name()))
                && key.contains(&format!("{}^1", p2.name()))
                && key.contains(&format!("{}^3", p3.name()))
                && key.contains(&format!("{}^5", p4.name()))
                && key.contains(&format!("{}^4", p5.name())),
            "hash_key() must contain perturbation orders"
        );

        let display = format!("{}", chain);
        assert!(
            display.contains(&format!("{}^2", p1.name()))
                && display.contains(&format!("{}^1", p2.name()))
                && display.contains(&format!("{}^3", p3.name()))
                && display.contains(&format!("{}^5", p4.name()))
                && display.contains(&format!("{}^4", p5.name())),
            "Display should contain correct names and orders"
        );
    }

    #[test]
    fn test_serialization() {
        let chain = make_pert_multichain(2u32, 8u32, 0u32, 10u32);
        let serialized = serde_json::to_string(&chain).unwrap();
        let deserialized: PertMultichain = serde_json::from_str(&serialized).unwrap();

        assert_eq!(chain, deserialized);
    }
}
