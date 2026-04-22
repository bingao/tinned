use std::collections::BTreeMap;
use std::fmt;
use std::sync::Arc;

use crate::perturbations::{PertMultichain, Perturbation};

pub trait PertSequence {
    fn as_vec(&self) -> Vec<Arc<Perturbation>>;
    fn as_map(&self) -> BTreeMap<Arc<Perturbation>, u32>;

    fn fmt_sequence(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, p) in self.as_vec().iter().enumerate() {
            if i > 0 {
                f.write_str(";")?;
            }
            write!(f, "{p}")?;
        }
        Ok(())
    }
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

impl fmt::Display for dyn PertSequence + '_ {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt_sequence(f)
    }
}
