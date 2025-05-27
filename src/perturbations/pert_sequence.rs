use std::collections::BTreeMap;
use std::sync::Arc;

use crate::perturbations::{PertMultichain, Perturbation};

pub trait PertSequence {
    fn as_vec(&self) -> Vec<Arc<Perturbation>>;
    fn as_map(&self) -> BTreeMap<Arc<Perturbation>, u32>;
}

impl PertSequence for &[Arc<Perturbation>] {
    #[inline]
    fn as_vec(&self) -> Vec<Arc<Perturbation>> {
        self.to_vec()
    }

    #[inline]
    fn as_map(&self) -> BTreeMap<Arc<Perturbation>, u32> {
        let mut map = BTreeMap::new();
        for p in self.iter() {
            *map.entry(p.clone()).or_insert(0) += 1;
        }

        map
    }
}

impl PertSequence for Vec<Arc<Perturbation>> {
    #[inline]
    fn as_vec(&self) -> Vec<Arc<Perturbation>> {
        self.clone()
    }

    #[inline]
    fn as_map(&self) -> BTreeMap<Arc<Perturbation>, u32> {
        let mut map = BTreeMap::new();
        for p in self {
            *map.entry(p.clone()).or_insert(0) += 1;
        }

        map
    }
}

impl PertSequence for PertMultichain {
    #[inline]
    fn as_vec(&self) -> Vec<Arc<Perturbation>> {
        self.to_vec()
    }

    #[inline]
    fn as_map(&self) -> BTreeMap<Arc<Perturbation>, u32> {
        self.get_map_clone()
    }
}
