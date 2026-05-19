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

    /// Creates a new perturbation multichain from a given slice of `Arc<Perturbation>`.
    #[inline]
    pub fn from_slice(slice: &[Arc<Perturbation>]) -> Self {
        let mut map = BTreeMap::new();
        for p in slice {
            *map.entry(p.clone()).or_insert(0) += 1;
        }
        PertMultichain(Arc::new(Mutex::new(map)))
    }

    /// Creates a new perturbation multichain by deeply cloning the underlying
    /// `BTreeMap` and inserting a given perturbation.
    #[inline]
    pub fn with_added_perturbation(&self, p: Arc<Perturbation>) -> Self {
        let mut map = self.get_map_clone();
        *map.entry(p).or_insert(0) += 1;
        Self::from_map(map)
    }

    /// Inserts a perturbation into the multichain, or increases the order by 1
    /// if the perturbation already exists in the multichain.
    #[inline]
    pub fn insert(&mut self, p: Arc<Perturbation>) {
        let mut map = self.0.lock().unwrap();
        let entry = map.entry(p).or_insert(0);
        *entry += 1;
    }

    //    /// Merges another perturbation multichain into this one in place.
    //    ///
    //    /// For each perturbation in `other`, its differentiation order is added to the
    //    /// corresponding order in `self`. If the perturbation does not already exist in
    //    /// `self`, it is inserted with the order from `other`.
    //    ///
    //    /// This method modifies `self`.
    //    ///
    //    /// # Panics
    //    ///
    //    /// Panics if the internal mutex of either `self` or `other` is poisoned.
    //    pub fn merge_in_place(&mut self, other: &PertMultichain) {
    //        let mut self_map = self.0.lock().unwrap();
    //
    //        let other_map = other.0.lock().unwrap();
    //
    //        for (perturbation, order) in other_map.iter() {
    //            let entry = self_map.entry(Arc::clone(perturbation)).or_insert(0);
    //
    //            *entry += order;
    //        }
    //    }
    //
    //    /// Returns a new perturbation multichain formed by merging `self` and `other`.
    //    ///
    //    /// For each perturbation appearing in either multichain, the resulting
    //    /// differentiation order is the sum of its orders in `self` and `other`.
    //    /// Perturbations that appear in only one multichain are copied unchanged into
    //    /// the result.
    //    ///
    //    /// This method does not modify either `self` or `other`.
    //    ///
    //    /// # Panics
    //    ///
    //    /// Panics if the internal mutex of either `self` or `other` is poisoned.
    //    #[inline]
    //    pub fn merge(&self, other: &PertMultichain) -> PertMultichain {
    //        let mut result_map = self.get_map_clone();
    //
    //        for (perturbation, order) in other.get_map_clone() {
    //            let entry = result_map.entry(perturbation).or_insert(0);
    //
    //            *entry += order;
    //        }
    //
    //        Self::from_map(result_map)
    //    }

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

    /// Returns the sum of all perturbation orders.
    #[inline]
    pub fn total_order(&self) -> u32 {
        self.0.lock().unwrap().values().sum()
    }

    /// Returns all perturbations in the multichain and meanwhile preserves the order.
    #[inline]
    pub fn keys(&self) -> Vec<Arc<Perturbation>> {
        let map = self.0.lock().unwrap();
        map.keys().cloned().collect()
    }

    /// Returns a cloned internal map of the perturbation multichain.
    #[inline]
    pub fn get_map_clone(&self) -> BTreeMap<Arc<Perturbation>, u32> {
        self.0.lock().unwrap().clone()
    }

    /// Returns a cloned vector of the perturbation multichain.
    #[inline]
    pub fn to_vec(&self) -> Vec<Arc<Perturbation>> {
        let map = self.0.lock().unwrap();
        let mut result = Vec::new();

        for (pert, &count) in map.iter() {
            if count > 0 {
                result.extend(std::iter::repeat(pert.clone()).take(count as usize));
            }
        }

        result
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

    /// Similar to the function `is_subchain` but takes the `subchain` in a
    /// vector of `Arc<Perturbation>`.
    #[inline]
    pub fn is_subchain_vec(&self, subchain: &[Arc<Perturbation>]) -> bool {
        let submap: BTreeMap<_, u32> = {
            let mut pert_map = BTreeMap::new();
            for pert in subchain {
                *pert_map.entry(pert.clone()).or_insert(0) += 1;
            }
            pert_map
        };

        let map = self.0.lock().unwrap();

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

    /// Similar to the function `is_superchain` but takes the `superchain` in a
    /// vector of `Arc<Perturbation>`.
    #[inline]
    pub fn is_superchain_vec(&self, superchain: &[Arc<Perturbation>]) -> bool {
        let supermap: BTreeMap<_, u32> = {
            let mut pert_map = BTreeMap::new();
            for pert in superchain {
                *pert_map.entry(pert.clone()).or_insert(0) += 1;
            }
            pert_map
        };

        let map = self.0.lock().unwrap();

        map.iter().all(|(p, &order)| supermap.get(p).copied().unwrap_or(0) >= order)
    }

    /// Checks if two `PertMultichains` share any common `Arc<Perturbation>`
    /// keys (i.e., their keys intersect).
    #[inline]
    pub fn has_overlap(&self, other: &PertMultichain) -> bool {
        let self_map = self.0.lock().unwrap().clone();
        let other_map = other.0.lock().unwrap().clone();

        self_map.keys().any(|k| other_map.contains_key(k))
    }

    /// Returns the "difference" between `self` and `other` as a
    /// `Vec<Arc<Perturbation>>`, where:
    ///
    /// (1) Each perturbation in `self` but not in `other` is added count times.
    /// (2) Each perturbation in both, but where `self` has a higher count, is
    ///     added (`self_count` - `other_count`) times.
    /// (3) Any perturbation only in `other` is ignored.
    #[inline]
    pub fn complement(&self, other: &PertMultichain) -> Vec<Arc<Perturbation>> {
        let self_map = self.0.lock().unwrap().clone();
        let other_map = other.0.lock().unwrap().clone();

        let mut result = Vec::new();

        for (pert, &count_self) in self_map.iter() {
            let count_other = other_map.get(pert).copied().unwrap_or(0);
            if count_self > count_other {
                result.extend(
                    std::iter::repeat(pert.clone()).take((count_self - count_other) as usize),
                );
            }
        }

        result
    }

    /// Generates a compact string suitable for hashing a perturbation multichain.
    /// Derived Ord for Perturbation makes sure the hash key of PertMultichain is stable.
    #[inline]
    pub fn hash_key(&self) -> String {
        let map = self.0.lock().unwrap();
        let mut parts = Vec::with_capacity(map.len());
        for (pert, order) in map.iter() {
            parts.push(format!("{}^{}", pert.hash_key(), order));
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
    pub fn make_pert_vec(len_name: u32, val_range: u32) -> Vec<Arc<Perturbation>> {
        vec![
            make_perturbation_i64(len_name, val_range),
            make_perturbation_f64(len_name, val_range),
            make_perturbation_complex(len_name, val_range),
            make_perturbation_rational(len_name, val_range),
            make_perturbation_symbol(len_name, val_range),
        ]
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

        chain.insert(p1.clone());
        assert_eq!(chain.get_order(&p1), 1);

        chain.insert(p1.clone());
        assert_eq!(chain.get_order(&p1), 2);

        let mut keys: Vec<Arc<Perturbation>> = chain.keys();

        assert_eq!(keys, vec![p1.clone()]);

        let p2 = make_perturbation_f64(2u32, 10u32);
        let p3 = make_perturbation_complex(2u32, 10u32);
        let p4 = make_perturbation_rational(2u32, 10u32);
        let p5 = make_perturbation_symbol(2u32, 4u32);

        chain.insert(p2.clone());
        chain.insert(p3.clone());
        chain.insert(p4.clone());
        chain.insert(p5.clone());

        keys = chain.keys();

        let mut expected_keys = vec![p1, p2, p3, p4, p5];
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

        chain.insert(p1.clone());
        chain.insert(p1.clone());
        chain.insert(p1.clone());
        chain.insert(p2.clone());
        chain.insert(p2.clone());
        chain.insert(p3.clone());
        chain.insert(p4.clone());

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
            c1.insert(p);
        }

        assert_eq!(c1, c2);

        assert!(c1.is_superchain(&c2));
        assert!(c1.is_subchain(&c2));
        assert!(c2.is_superchain(&c1));
        assert!(c2.is_subchain(&c1));

        c1.insert(make_perturbation_i64(2u32, 10u32));
        c2.insert(make_perturbation_f64(2u32, 10u32));

        assert!(!c1.is_superchain(&c2));
        assert!(!c1.is_subchain(&c2));
        assert!(!c2.is_superchain(&c1));
        assert!(!c2.is_subchain(&c1));

        let c3 = c1.keys();

        assert!(c1.is_subchain_vec(&c3));
    }

    #[test]
    fn test_with_added_perturbation() {
        let mut c1 = make_pert_multichain(2u32, 8u32, 0u32, 10u32);
        let p = make_perturbation_symbol(2u32, 4u32);
        let c2 = c1.with_added_perturbation(p.clone());

        assert!(c1.is_superchain(&c2));

        c1.insert(p);

        assert_eq!(c1, c2);
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
            key.contains(&format!("{}^2", p1.hash_key()))
                && key.contains(&format!("{}^1", p2.hash_key()))
                && key.contains(&format!("{}^3", p3.hash_key()))
                && key.contains(&format!("{}^5", p4.hash_key()))
                && key.contains(&format!("{}^4", p5.hash_key())),
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
