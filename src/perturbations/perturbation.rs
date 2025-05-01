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
        crate::internal::intern_pert(Arc::new(Self {
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
pub mod test_utils {
    use super::*;
    use crate::expressions::number::test_utils::*;
    use crate::expressions::symbol::test_utils::*;

    #[inline]
    pub fn make_perturbation_i64(len_name: u32, val_range: u32) -> Arc<Perturbation> {
        let name: String = if len_name == 0 {
            "a".to_string()
        } else {
            random_alphanumeric(len_name)
        };
        Perturbation::new(name, make_number_i64(val_range))
    }

    #[inline]
    pub fn make_perturbation_f64(len_name: u32, val_range: u32) -> Arc<Perturbation> {
        let name: String = if len_name == 0 {
            "a".to_string()
        } else {
            random_alphanumeric(len_name)
        };
        Perturbation::new(name, make_number_f64(val_range))
    }

    #[inline]
    pub fn make_perturbation_complex(len_name: u32, val_range: u32) -> Arc<Perturbation> {
        let name: String = if len_name == 0 {
            "a".to_string()
        } else {
            random_alphanumeric(len_name)
        };
        Perturbation::new(name, make_number_complex(val_range))
    }

    #[inline]
    pub fn make_perturbation_rational(len_name: u32, val_range: u32) -> Arc<Perturbation> {
        let name: String = if len_name == 0 {
            "a".to_string()
        } else {
            random_alphanumeric(len_name)
        };
        Perturbation::new(name, make_number_rational(val_range))
    }

    #[inline]
    pub fn make_perturbation_symbol(len_name: u32, len_freq: u32) -> Arc<Perturbation> {
        let name: String = if len_name == 0 {
            "a".to_string()
        } else {
            random_alphanumeric(len_name)
        };
        Perturbation::new(name, make_symbol(len_freq))
    }
}

#[cfg(test)]
mod tests {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    use super::test_utils::*;
    use super::*;
    use crate::expressions::Number;

    test_struct_safety!(Perturbation);

    test_thread_interning!({ make_perturbation_complex(0u32, 0u32) });

    #[test]
    fn test_struct() {
        let f1: Arc<dyn Expr> = Number::Integer(1).into();
        let p1 = Perturbation::new("alpha", f1.clone());

        assert_eq!(p1.name(), "alpha");
        assert_eq!(p1.frequency(), &f1);
        assert!(p1.to_string().contains("alpha"));

        let f2: Arc<dyn Expr> = Number::Real(2.0).into();
        let p2 = Perturbation::new("beta", f2);
        let p3 = Perturbation::new("alpha", f1);

        assert_ne!(p1, p2);
        assert_eq!(p1, p3);
    }

    #[test]
    fn test_serialization() {
        let mut p = make_perturbation_i64(2u32, 10u32);
        let mut serialized = serde_json::to_string(&p).unwrap();
        let mut deserialized: Arc<Perturbation> = serde_json::from_str(&serialized).unwrap();
        assert_eq!(p, deserialized);

        p = make_perturbation_f64(2u32, 10u32);
        serialized = serde_json::to_string(&p).unwrap();
        deserialized = serde_json::from_str(&serialized).unwrap();
        assert_eq!(p, deserialized);

        p = make_perturbation_complex(2u32, 10u32);
        serialized = serde_json::to_string(&p).unwrap();
        deserialized = serde_json::from_str(&serialized).unwrap();
        assert_eq!(p, deserialized);

        p = make_perturbation_rational(2u32, 10u32);
        serialized = serde_json::to_string(&p).unwrap();
        deserialized = serde_json::from_str(&serialized).unwrap();
        assert_eq!(p, deserialized);

        p = make_perturbation_symbol(2u32, 4u32);
        serialized = serde_json::to_string(&p).unwrap();
        deserialized = serde_json::from_str(&serialized).unwrap();
        assert_eq!(p, deserialized);
    }

    #[test]
    fn test_interning() {
        let p1 = make_perturbation_complex(2u32, 10u32);
        let p2 = Perturbation::new(p1.name(), p1.frequency().clone());

        assert!(Arc::ptr_eq(&p1, &p2));
    }

    #[test]
    fn test_order_and_hash() {
        let p1 = make_perturbation_complex(1u32, 10u32);
        let p2 = make_perturbation_symbol(2u32, 10u32);

        let p3 = Perturbation::new(p1.name(), p2.frequency().clone());

        assert!(p1 < p3 || p1 > p3);

        if p1.name().to_string() <= p2.name().to_string() {
            assert!(p1 <= p2);
        } else {
            assert!(p1 > p2);
        }

        let mut hasher1 = DefaultHasher::new();
        p1.hash(&mut hasher1);
        let mut hasher2 = DefaultHasher::new();
        p1.hash(&mut hasher2);
        assert_eq!(hasher1.finish(), hasher2.finish());
    }
}
