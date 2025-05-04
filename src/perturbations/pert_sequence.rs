use std::sync::Arc;

use crate::perturbations::{PertMultichain, Perturbation};

pub trait PertSequence {
    fn ordered_perturbations(&self) -> Vec<Arc<Perturbation>>;
}

impl PertSequence for Vec<Arc<Perturbation>> {
    fn ordered_perturbations(&self) -> Vec<Arc<Perturbation>> {
        self.clone()
    }
}

impl PertSequence for PertMultichain {
    fn ordered_perturbations(&self) -> Vec<Arc<Perturbation>> {
        self.to_vec()
    }
}
