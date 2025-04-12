use std::sync::Arc;

use crate::core::Expr;

// Perturbations are compared first by name (String), then by frequency
// (Arc<dyn Expr>).
//
// So two perturbations with the same name but different frequencies will not
// be equal and will have a deterministic ordering.
#[derive(Clone, Debug, PartialOrd, Ord, Hash, serde::Serialize, serde::Deserialize)]
pub struct Perturbation {
    name: String,
    frequency: Arc<dyn Expr>,
}

impl Perturbation {
    #[inline]
    pub fn new(name: impl Into<String>, frequency: Arc<dyn Expr>) -> Arc<Self> {
        crate::utils::intern_pert(Arc::new(Self {
            name: name.into(),
            frequency,
        }))
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[inline]
    pub fn frequency(&self) -> &Arc<dyn Expr> {
        &self.frequency
    }

    // A canonical string used for hashing, ordering, interning, etc.
    #[inline]
    pub fn hash_key(&self) -> String {
        format!("Perturbation({},{})", self.name, self.frequency.hash_key())
    }
}

impl PartialEq for Perturbation {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name && &self.frequency == &other.frequency
    }
}

impl Eq for Perturbation {}

impl std::fmt::Display for Perturbation {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}({})", self.name, self.frequency)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    use serde_json;

    use super::*;
    use crate::expressions::{Number, Symbol};

    #[test]
    fn test_struct() {
        let freq: Arc<dyn Expr> = Number::Integer(1).into();
        let pert = Perturbation::new("alpha", freq.clone());

        assert_eq!(pert.name(), "alpha");
        assert_eq!(pert.frequency(), &freq);
        assert!(pert.to_string().contains("alpha"));
    }

    #[test]
    fn test_serialization() {
        let freq: Arc<dyn Expr> = Symbol::new("omega");
        let pert = Perturbation::new("beta", freq.clone());

        let serialized = serde_json::to_string(&pert).unwrap();
        let deserialized: Arc<Perturbation> = serde_json::from_str(&serialized).unwrap();

        assert_eq!(pert, deserialized);
    }

    #[test]
    fn test_interning() {
        let f: Arc<dyn Expr> = Number::Integer(5).into();
        let p1 = Perturbation::new("gamma", f.clone());
        let p2 = Perturbation::new("gamma", f);

        assert!(Arc::ptr_eq(&p1, &p2));
    }

    #[test]
    fn test_order_and_hash() {
        let a: Arc<dyn Expr> = Symbol::new("x");
        let b: Arc<dyn Expr> = Symbol::new("y");

        let p1 = Perturbation::new("a", a.clone());
        let p2 = Perturbation::new("a", b.clone());
        let p3 = Perturbation::new("b", a.clone());

        // Test Ord
        assert!(p1 < p2 || p1 > p2); // Different freq
        assert!(p1 < p3); // Different name

        // Test Hash
        let mut hasher1 = DefaultHasher::new();
        p1.hash(&mut hasher1);
        let mut hasher2 = DefaultHasher::new();
        p1.hash(&mut hasher2);
        assert_eq!(hasher1.finish(), hasher2.finish());
    }
}
