use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};

// Exponential adjoint map (or conjugation operation in Lie algebra):
// exp(ad_{X})(Y) = exp(X)*Y*exp(-X) (`left_action` is true), or
// exp(ad_{-X})(Y) = exp(-X)*Y*exp(X) (`left_action` is false).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct ExpAdjointMap {
    generator: Arc<dyn Expr>,
    generator_commutative: bool,
    target: Arc<dyn Expr>,
    left_action: bool,
}

impl ExpAdjointMap {

}

// Generate all partitions of {1, ..., n} into exactly k non-empty subsets
fn generate_partitions(n: usize, k: usize) -> Vec<Vec<Vec<usize>>> {
    if n < k || k == 0 {
        return vec![]; // No valid partitions
    }

    if n == k {
        return vec![ (1..=n).map(|i| vec![i]).collect() ]; // Each element in its own group
    }

    if k == 1 {
        return vec![ vec![(1..=n).collect()] ]; // All elements in one group
    }

    let mut new_partitions = Vec::new();

    // Case 1: Use partitions from S(n-1, k), placing n into existing groups
    for partition in generate_partitions(n - 1, k) {
        for i in 0..partition.len() {
            let mut new_partition = partition.clone();
            new_partition[i].push(n);
            new_partitions.push(new_partition);
        }
    }

    // Case 2: Use partitions from S(n-1, k-1), creating a new singleton {n}
    for partition in generate_partitions(n - 1, k - 1) {
        let mut new_partition = partition.clone();
        new_partition.push(vec![n]);
        new_partitions.push(new_partition);
    }

    new_partitions
}

impl ExprInternal for ExpAdjointMap {
    impl_expr_internal_methods!(ExpAdjointMap);

    #[inline]
    fn hash_key(&self) -> String {
    }

    #[inline]
    fn match_for_find_all(&self, other: &Arc<dyn Expr>) -> bool {
    }
}

#[typetag::serde]
impl Expr for ExpAdjointMap {
    impl_expr_common_methods!(false);

    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {

    }
}

impl PartialEq for ExpAdjointMap {
    fn eq(&self, other: &Self) -> bool {
    }
}

impl Eq for ExpAdjointMap {}

impl std::fmt::Display for ExpAdjointMap {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    test_struct_safety!(ExpAdjointMap);
}
