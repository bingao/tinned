use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

use serde::{Serialize, Deserialize, Serializer, Deserializer};

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
        let map = self.0.lock().unwrap();
        let submap = subchain.0.lock().unwrap();

        map.iter().all(|(p, &order)| submap.get(p).copied().unwrap_or(0) <= order)
    }

    /// Returns true if `chain` is a sub-multichain of `superchain`. That is,
    /// all perturbations in `chain` appear in `superchain` with at least the
    /// same order.
    #[inline]
    pub fn is_superchain(&self, superchain: &PertMultichain) -> bool {
        let map = self.0.lock().unwrap();
        let supermap = superchain.0.lock().unwrap();

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
        let self_map = self.0.lock().unwrap();
        let other_map = other.0.lock().unwrap();
        *self_map == *other_map
    }
}

impl Eq for PertMultichain {}

impl std::fmt::Display for PertMultichain {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let map = self.0.lock().unwrap();
        let parts: Vec<String> =
            map.iter().map(|(pert, order)| format!("{}^{}", pert.name(), order)).collect();
        f.write_str(&parts.join(", "))
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
        let entries: Vec<PertEntry> = map.iter()
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
mod tests {
    use super::*;
    use crate::expressions::{Number, Symbol};
    use crate::perturbations::Perturbation;

    test_struct_safety!(PertMultichain);

    #[test]
    fn test_insert_and_get_order() {
        let mut chain = PertMultichain::new();
        let p = Perturbation::new("alpha", Number::from_i64(0));

        assert_eq!(chain.get_order(&p), 0);

        chain.insert(&p);
        assert_eq!(chain.get_order(&p), 1);

        chain.insert(&p);
        assert_eq!(chain.get_order(&p), 2);
    }

    #[test]
    fn test_get_map_clone() {
        let mut chain = PertMultichain::new();
        let p1 = Perturbation::new("p1", Symbol::new("x"));
        let p2 = Perturbation::new("p2", Number::from_i64(42));

        chain.insert(&p1);
        chain.insert(&p1);
        chain.insert(&p2);

        let map = chain.get_map_clone();

        assert_eq!(map.get(&p1), Some(&2));
        assert_eq!(map.get(&p2), Some(&1));
    }

    #[test]
    fn test_is_subchain_and_superchain() {
        let mut c1 = PertMultichain::new();
        let mut c2 = PertMultichain::new();

        let p1 = Perturbation::new("p1", Symbol::new("x"));
        let p2 = Perturbation::new("p2", Number::from_i64(42));

        c1.insert(&p1);
        c1.insert(&p1);
        c1.insert(&p2);

        c2.insert(&p1);

        assert!(!c1.is_superchain(&c2));
        assert!(c1.is_subchain(&c2));
        assert!(c2.is_superchain(&c1));
        assert!(!c2.is_subchain(&c1));

        c2.insert(&p2);

        assert!(!c1.is_superchain(&c2));
        assert!(c1.is_subchain(&c2));
        assert!(c2.is_superchain(&c1));
        assert!(!c2.is_subchain(&c1));

        c2.insert(&p2);

        assert!(!c1.is_superchain(&c2));
        assert!(!c1.is_subchain(&c2));
        assert!(!c2.is_superchain(&c1));
        assert!(!c2.is_subchain(&c1));
    }

    #[test]
    fn test_hash_key_and_display() {
        let mut chain = PertMultichain::new();
        let p1 = Perturbation::new("alpha", Symbol::new("x"));
        let p2 = Perturbation::new("beta", Number::from_i64(5));

        chain.insert(&p1);
        chain.insert(&p2);
        chain.insert(&p1);

        let key = chain.hash_key();
        assert!(
            key.contains("alpha^2") && key.contains("beta^1"),
            "hash_key() must contain perturbation orders"
        );

        let display = format!("{}", chain);
        assert!(
            display.contains("alpha^2") && display.contains("beta^1"),
            "Display should contain correct names and orders"
        );
    }

    #[test]
    fn test_equality() {
        let p1 = Perturbation::new("p1", Symbol::new("x"));
        let p2 = Perturbation::new("p2", Number::from_i64(42));

        let mut c1 = PertMultichain::new();
        let mut c2 = PertMultichain::new();

        c1.insert(&p1);
        c1.insert(&p1);
        c1.insert(&p2);

        c2.insert(&p1);
        c2.insert(&p1);
        c2.insert(&p2);

        assert_eq!(c1, c2);
    }

    #[test]
    fn test_serialization() {
        use serde_json;

        let p1 = Perturbation::new("p1", Symbol::new("x"));
        let p2 = Perturbation::new("p2", Number::from_i64(42));

        let mut c = PertMultichain::new();
        c.insert(&p1);
        c.insert(&p1);
        c.insert(&p2);

        let serialized = serde_json::to_string(&c).unwrap();
        let deserialized: PertMultichain = serde_json::from_str(&serialized).unwrap();

        assert_eq!(c, deserialized);
    }
}
