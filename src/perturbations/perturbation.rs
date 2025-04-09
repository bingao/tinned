use std::cmp::Ordering;
use std::sync::Arc;

use crate::core::Expr;

#[derive(
    Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, serde::Serialize, serde::Deserialize,
)]
pub struct Perturbation {
    name: String,
    frequency: Arc<dyn Expr>,
}

impl Perturbation {
    #[inline]
    pub fn new(name: impl Into<String>, frequency: Arc<dyn Expr>) -> Arc<Self> {
        //crate::utils::intern(Arc::new(Self { name: name.into(), frequency }))
        Arc::new(Self { name: name.into(), frequency })
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[inline]
    pub fn frequency(&self) -> &Arc<dyn Expr> {
        &self.frequency
    }

    // A canonical string used for hashing and ordering
    #[inline]
    pub fn hash_key(&self) -> String {
        format!("Perturbation({},{})", self.name, self.frequency.hash_key())
    }
}

impl PartialOrd for Perturbation {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.hash_key().cmp(&other.hash_key()))
    }
}

impl Ord for Perturbation {
    fn cmp(&self, other: &Self) -> Ordering {
        self.hash_key().cmp(&other.hash_key())
    }
}

impl std::fmt::Display for Perturbation {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}({})", self.name, self.frequency)
    }
}
