use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{MatrixAdd, ZeroOperator};
use crate::internal::{
    intern_expr, multi_expression_format, multi_expression_hash, sort_expressions_grouped_by,
};
use crate::perturbations::Perturbation;
use crate::public::{
    NumberTolerance, downcast_from_arc, downcast_from_ref, generic_expression_error, is_zero_expr,
};

// Either adjoint map (or adjoint action, adjoint representation)
// [xn, [..., [x1, [x0, y]]...]], or its "right" acting version
// [[...[[y, x0], x1], ...], xn]. Users decide which version they are using.
// The default is acting from the left side.
//
// `adjoint_chain` holds x0, x1, ..., xn, and y is stored in `target`. The
// field `chain_commutative` indicates if [xi, xj] = 0 for all 0 <= i, j <= n.
// That is, if `chain_commutative` is `true`, then we can sort `adjoint_chain`
// in some way without changing the result of adjoint map.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct AdjointMap {
    adjoint_chain: Vec<Arc<dyn Expr>>,
    chain_commutative: bool,
    target: Arc<dyn Expr>,
    left_action: bool,
}

impl AdjointMap {
    pub fn new(
        adjoint_chain: Vec<Arc<dyn Expr>>,
        chain_commutative: bool,
        target: Arc<dyn Expr>,
        left_action: Option<bool>,
    ) -> Arc<dyn Expr> {
        let left_action = left_action.unwrap_or(true);

        return if chain_commutative {
            // Sort the chain according to `total_order()`
            let sorted = sort_expressions_grouped_by(&adjoint_chain, |e| e.total_order());
            intern_expr(Arc::new(Self {
                adjoint_chain: sorted,
                chain_commutative,
                target,
                left_action,
            }))
        } else {
            intern_expr(Arc::new(Self {
                adjoint_chain,
                chain_commutative,
                target,
                left_action,
            }))
        };
    }

    #[inline]
    pub fn adjoint_chain(&self) -> &[Arc<dyn Expr>] {
        &self.adjoint_chain
    }

    #[inline]
    pub fn chain_commutative(&self) -> bool {
        self.chain_commutative
    }

    #[inline]
    pub fn target(&self) -> &Arc<dyn Expr> {
        &self.target
    }

    #[inline]
    pub fn left_action(&self) -> bool {
        self.left_action
    }
}

impl ExprInternal for AdjointMap {
    impl_expr_internal_methods!(AdjointMap);

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "AdjointMap([{}]; {}; {}; {})",
            multi_expression_hash(&self.adjoint_chain, ";"),
            self.chain_commutative,
            self.target.hash_key(),
            self.left_action
        )
    }

    #[inline]
    fn match_for_find_all(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<AdjointMap>(other) {
            let len = self.adjoint_chain.len();

            // We find adjoint maps with `left_action` either true or false
            if !self.target.match_for_find_all(&op.target) || len != op.adjoint_chain.len() {
                return false;
            }

            match (self.chain_commutative, op.chain_commutative) {
                (false, false) | (true, true) => self
                    .adjoint_chain
                    .iter()
                    .zip(&op.adjoint_chain)
                    .all(|(a, b)| a.match_for_find_all(b)),
                (true, false) => {
                    let sorted =
                        sort_expressions_grouped_by(&op.adjoint_chain, |e| e.total_order());
                    self.adjoint_chain.iter().zip(&sorted).all(|(a, b)| a.match_for_find_all(b))
                },
                (false, true) => {
                    let sorted =
                        sort_expressions_grouped_by(&self.adjoint_chain, |e| e.total_order());
                    sorted.iter().zip(&op.adjoint_chain).all(|(a, b)| a.match_for_find_all(b))
                },
            }
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for AdjointMap {
    impl_expr_common_methods!(false);

    #[inline]
    fn clean_temporum(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        impl_adjoint_map_operation!(
            self,
            |x: &Arc<dyn Expr>| x.clean_temporum(freq_tol.clone()),
            "AdjointMap::clean_temporum() failed"
        )
    }

    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        // Precompute the derivative of each x and store it
        let with_context = |f: &Arc<dyn Expr>| {
            f.differentiate(s).map_err(|e| {
                generic_expression_error(
                    "AdjointMap::differentiate() failed for chain",
                    self,
                    Some(Box::new(e)),
                )
            })
        };

        let diff_adj_chain: Vec<Arc<dyn Expr>> =
            self.adjoint_chain.iter().map(with_context).collect::<Result<_, _>>()?;

        let mut results = Vec::with_capacity(diff_adj_chain.len() + 1);

        for (i, diff) in diff_adj_chain.iter().enumerate() {
            // Skip derivative = 0 to avoid 0 * [...] = 0
            if is_zero_expr(diff, None) {
                continue;
            }

            let mut new_adj_chain = self.adjoint_chain.clone();
            // For each x, replace it with its derivative while keeping others
            // intact
            new_adj_chain[i] = diff.clone();

            results.push(Self::new(
                new_adj_chain,
                self.chain_commutative,
                self.target.clone(),
                Some(self.left_action),
            ));
        }

        let diff_target = self.target.differentiate(s).map_err(|e| {
            generic_expression_error(
                "AdjointMap::differentiate() failed for target",
                self,
                Some(Box::new(e)),
            )
        })?;

        if !is_zero_expr(&diff_target, None) {
            results.push(Self::new(
                self.adjoint_chain.clone(),
                self.chain_commutative,
                diff_target,
                Some(self.left_action),
            ));
        }

        MatrixAdd::new(results)
    }

    #[inline]
    fn eliminate(
        &self,
        parameter: &Arc<dyn Expr>,
        perturbations: &[Arc<Perturbation>],
        min_order: u32,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        impl_adjoint_map_operation!(
            self,
            |x: &Arc<dyn Expr>| x.eliminate(parameter, perturbations, min_order),
            "AdjointMap::eliminate() failed"
        )
    }

    #[inline]
    fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
        if self.adjoint_chain.iter().any(|x| x.exist_any(set)) {
            return true;
        }

        set.iter().any(|expr| self.eq_expr(expr.as_ref())) || self.target.exist_any(set)
    }

    #[inline]
    fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
        if self.match_for_find_all(s) {
            return BTreeMap::from([(self.total_order(), HashSet::from([self.clone_expr()]))]);
        }

        let mut result = self.target.find_all(s);

        for x in &self.adjoint_chain {
            for (order, subset) in x.find_all(s) {
                result.entry(order).or_default().extend(subset);
            }
        }

        result
    }

    #[inline]
    fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
            return Ok(ZeroOperator::new());
        }

        impl_adjoint_map_operation!(
            self,
            |x: &Arc<dyn Expr>| x.remove(set),
            "AdjointMap::remove() failed"
        )
    }

    #[inline]
    fn replace(
        &self,
        map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if let Some((_, value)) = map.iter().find(|(key, _)| self.eq_expr(key.as_ref())) {
            return Ok(value.clone());
        }

        impl_adjoint_map_operation!(
            self,
            |x: &Arc<dyn Expr>| x.replace(map),
            "AdjointMap::replace() failed"
        )
    }

    #[inline]
    fn replace_all(
        &self,
        map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if let Some((_, value)) = map.iter().find(|(key, _)| self.match_for_replace_all(key)) {
            return Ok(value.clone());
        }

        impl_adjoint_map_operation!(
            self,
            |x: &Arc<dyn Expr>| x.replace_all(map),
            "AdjointMap::replace_all() failed"
        )
    }
}

impl PartialEq for AdjointMap {
    fn eq(&self, other: &Self) -> bool {
        let len = self.adjoint_chain.len();

        if &self.target != &other.target
            || (self.left_action != other.left_action && len % 2 == 1)
            || len != other.adjoint_chain.len()
        {
            return false;
        }

        match (self.chain_commutative, other.chain_commutative) {
            (false, false) | (true, true) => self.adjoint_chain == other.adjoint_chain,
            (true, false) => {
                self.adjoint_chain
                    == sort_expressions_grouped_by(&other.adjoint_chain, |e| e.total_order())
            },
            (false, true) => {
                sort_expressions_grouped_by(&self.adjoint_chain, |e| e.total_order())
                    == other.adjoint_chain
            },
        }
    }
}

impl Eq for AdjointMap {}

impl std::fmt::Display for AdjointMap {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.left_action {
            write!(
                f,
                "({})[{},{}{}",
                self.chain_commutative,
                multi_expression_format(&self.adjoint_chain, ",["),
                self.target,
                "]".repeat(self.adjoint_chain.len())
            )
        } else {
            write!(
                f,
                "({}){}{},{}]",
                self.chain_commutative,
                "[".repeat(self.adjoint_chain.len()),
                self.target,
                multi_expression_format(&self.adjoint_chain, "],")
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    test_struct_safety!(AdjointMap);
}
