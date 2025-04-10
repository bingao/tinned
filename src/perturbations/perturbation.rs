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
        crate::utils::intern_pert(Arc::new(Self { name: name.into(), frequency }))
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
