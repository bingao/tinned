use std::collections::HashMap;
use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{
    Add, MatrixAdd, MatrixMul, Mul, Number, OneElecMatrix, TimeEvolution, ZeroOperator,
};
use crate::internal::intern_expr;
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{
    NumberTolerance, downcast_from_arc, generic_expression_error, is_expr_type, is_zero_expr,
    unreachable_error,
};

// Represents the operator (18), J. Comput. Chem. 2024; 45: 2136-2152.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct BasisTimeEvolution {
    at_zero_perturbations: bool,
    braket: Arc<dyn Expr>,
    dependencies: PertMultichain,
    derivative: PertMultichain,
}

impl BasisTimeEvolution {
    // Note: `dependencies` is the perturbation dependencies of Sb and Sk
    #[inline]
    pub fn builder(dependencies: PertMultichain) -> BasisTimeEvolutionBuilder {
        BasisTimeEvolutionBuilder {
            at_zero_perturbations: Some(false),
            braket: None,
            dependencies,
            derivative: None,
        }
    }

    #[inline]
    fn with_braket(
        &self,
        braket: Arc<dyn Expr>,
        at_zero_perturbations: Option<bool>,
    ) -> BasisTimeEvolutionBuilder {
        BasisTimeEvolutionBuilder {
            at_zero_perturbations,
            braket: Some(braket),
            dependencies: self.dependencies.clone(),
            derivative: Some(self.derivative.clone()),
        }
    }

    #[inline]
    fn with_derivative(
        &self,
        derivative: PertMultichain,
        braket: Arc<dyn Expr>,
    ) -> BasisTimeEvolutionBuilder {
        BasisTimeEvolutionBuilder {
            at_zero_perturbations: Some(self.at_zero_perturbations),
            braket: Some(braket),
            dependencies: self.dependencies.clone(),
            derivative: Some(derivative),
        }
    }

    #[inline]
    pub fn at_zero_perturbations(&self) -> bool {
        self.at_zero_perturbations
    }

    #[inline]
    pub fn braket(&self) -> &Arc<dyn Expr> {
        &self.braket
    }

    #[inline]
    pub fn dependencies(&self) -> &PertMultichain {
        &self.dependencies
    }

    #[inline]
    pub fn derivative(&self) -> &PertMultichain {
        &self.derivative
    }

    // Returns frequency factor, bra and ket of all terms in `braket`. The
    // frequency factor is computed by using Equation (62), J. Comput. Chem.
    // 2024; 45: 2136-2152.
    pub(crate) fn at_zero_strength(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Vec<(Arc<dyn Expr>, Arc<dyn Expr>, Arc<dyn Expr>)>, TinnedError> {
        let terms = if let Some(mat_add) = downcast_from_arc::<MatrixAdd>(&self.braket) {
            mat_add.terms()
        } else {
            std::slice::from_ref(&self.braket)
        };

        let mut result = Vec::new();

        for term in terms {
            let mat_mul = downcast_from_arc::<MatrixMul>(term)
                .ok_or_else(|| unreachable_error("Unexpected term inside braket", term, None))?;

            let factors = mat_mul.factors();
            if factors.len() != 2 {
                return Err(unreachable_error(
                    "Unexpected number of factors of term inside braket",
                    term,
                    None,
                ));
            }

            if self.at_zero_perturbations {
                result.push((
                    mat_mul.coefficient().clone(),
                    factors[0].clone(),
                    factors[1].clone(),
                ));
            } else {
                let bra = downcast_from_arc::<TimeEvolution>(&factors[0]).ok_or_else(|| {
                    unreachable_error(
                        "Unexpected first factor of term inside braket",
                        &factors[0],
                        None,
                    )
                })?;

                let ket = downcast_from_arc::<TimeEvolution>(&factors[1]).ok_or_else(|| {
                    unreachable_error(
                        "Unexpected second factor of term inside braket",
                        &factors[1],
                        None,
                    )
                })?;

                let frequency = Mul::new(vec![
                    Add::new(vec![bra.frequency()?, ket.frequency()?])?,
                    mat_mul.coefficient().clone(),
                    Number::one_half(),
                ])?;

                if is_zero_expr(&frequency, freq_tol.clone()) {
                    continue;
                }

                result.push((frequency, bra.argument().clone(), ket.argument().clone()));
            }
        }

        Ok(result)
    }
}

// Helper function to build `braket` with given dependencies, by following
// Equation (62), J. Comput. Chem. 2024; 45: 2136-2152.
fn build_braket(deps: &PertMultichain) -> Result<Arc<dyn Expr>, TinnedError> {
    let bra = OneElecMatrix::builder("Sb").dependencies(deps.clone()).build()?;
    let dt_bra = TimeEvolution::builder(bra).is_forward(true).build()?;

    let ket = OneElecMatrix::builder("Sk").dependencies(deps.clone()).build()?;
    let dt_ket = TimeEvolution::builder(ket).is_forward(false).build()?;

    MatrixMul::new(vec![dt_bra, dt_ket])
}

#[derive(Debug)]
pub struct BasisTimeEvolutionBuilder {
    at_zero_perturbations: Option<bool>,
    braket: Option<Arc<dyn Expr>>,
    dependencies: PertMultichain,
    derivative: Option<PertMultichain>,
}

impl BasisTimeEvolutionBuilder {
    //#[inline]
    //fn at_zero_perturbations(mut self, at_zero_perturbations: bool) -> Self {
    //    self.at_zero_perturbations = Some(at_zero_perturbations);
    //    self
    //}

    //#[inline]
    //fn braket(mut self, braket: Arc<dyn Expr>) -> Self {
    //    self.braket = Some(braket);
    //    self
    //}

    //#[inline]
    //fn derivative(mut self, derivative: PertMultichain) -> Self {
    //    self.derivative = Some(derivative);
    //    self
    //}

    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        let at_zero_perturbations = self.at_zero_perturbations.unwrap_or(false);
        let braket = self.braket.unwrap_or(build_braket(&self.dependencies)?);
        let derivative = self.derivative.unwrap_or(PertMultichain::new());

        Ok(intern_expr(Arc::new(BasisTimeEvolution {
            at_zero_perturbations,
            braket,
            dependencies: self.dependencies,
            derivative,
        })))
    }
}

impl ExprInternal for BasisTimeEvolution {
    impl_expr_internal_methods!(BasisTimeEvolution, true);

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "BasisTimeEvolution({}; [{}]; {}; [{}])",
            self.at_zero_perturbations,
            self.dependencies.hash_key(),
            self.braket.hash_key(),
            self.derivative.hash_key(),
        )
    }

    #[inline]
    fn total_order(&self) -> u32 {
        self.derivative.total_order()
    }

    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<BasisTimeEvolution>(other) {
            // We care only `dependencies`, regardless whether at zero strength
            // or not (specified by `at_zero_perturbations`, `braket` also changes)
            self.dependencies == op.dependencies && self.derivative.is_subchain(&op.derivative)
        } else {
            false
        }
    }

    #[inline]
    fn eq_by_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<BasisTimeEvolution>(other) {
            self.at_zero_perturbations == op.at_zero_perturbations
                && self.dependencies == op.dependencies
                && self.derivative.is_subchain(&op.derivative)
        } else {
            false
        }
    }

    #[inline]
    fn replace_expr_children(
        &self,
        _map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        _include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(self.clone_expr())
    }

    #[inline]
    fn retain_single(
        &self,
        s: &Arc<dyn Expr>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_self_single(s, include_derivatives) {
            Ok(self.clone_expr())
        } else {
            Ok(ZeroOperator::new())
        }
    }
}

#[typetag::serde]
impl Expr for BasisTimeEvolution {
    impl_nullary_expr_common_methods!(BasisTimeEvolution, false);

    #[inline]
    fn has_unperturbed_term(&self) -> bool {
        false
    }

    // `BasisTimeEvolution` will disappear if it is unperturbed or all
    // perturbations have zero frequency
    fn substitute_zero_perturbations(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.at_zero_perturbations {
            return Ok(self.clone_expr());
        } else if self.derivative.is_empty() {
            return Ok(ZeroOperator::new());
        }

        let triplets = self.at_zero_strength(freq_tol)?;

        if triplets.is_empty() {
            return Ok(ZeroOperator::new());
        }

        let mut terms: Vec<Arc<dyn Expr>> = Vec::with_capacity(triplets.len());
        for triplet in triplets {
            terms.push(MatrixMul::new(vec![triplet.0, triplet.1, triplet.2])?);
        }

        self.with_braket(MatrixAdd::new(terms)?, Some(true)).build()
    }

    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_braket = self.braket.differentiate(s).map_err(|e| {
            generic_expression_error(
                "BasisTimeEvolution::differentiate() failed",
                self,
                Some(Box::new(e)),
            )
        })?;

        if is_expr_type::<ZeroOperator>(&diff_braket) {
            return Ok(diff_braket);
        }

        let new_deriv = self.derivative.with_added_perturbation(s);

        self.with_derivative(new_deriv, diff_braket).build()
    }

    // `BasisTimeEvolution` is an undivided whole for methods `exist_any()`,
    // `find_superchains()`. So, we use the corresponding methods of the
    // pub trait `Expr`.
}

impl PartialEq for BasisTimeEvolution {
    fn eq(&self, other: &Self) -> bool {
        self.at_zero_perturbations == other.at_zero_perturbations
            && &self.braket == &other.braket
            && self.dependencies == other.dependencies
            && self.derivative == other.derivative
    }
}

impl Eq for BasisTimeEvolution {}

impl std::fmt::Display for BasisTimeEvolution {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.derivative.is_empty() {
            write!(f, "op(T)")
        } else {
            write!(f, "op(T)^{}", self.derivative)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perturbations::pert_multichain::test_utils::{
        make_pert_multichain, make_super_multichain,
    };
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::is_one_expr;

    test_struct_safety!(BasisTimeEvolution);

    test_thread_interning!(
        BasisTimeEvolution::builder(make_pert_multichain(0u32, 0u32, 1u32, 0u32)).build().unwrap()
    );

    #[test]
    fn test_impl_expr() {
        let deps = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
        let op1 = BasisTimeEvolution::builder(deps.clone()).build().unwrap();

        let op = downcast_from_arc::<BasisTimeEvolution>(&op1).unwrap();
        let braket = build_braket(&deps).unwrap();
        assert_eq!(
            op,
            &BasisTimeEvolution {
                at_zero_perturbations: false,
                braket: braket.clone(),
                dependencies: deps.clone(),
                derivative: PertMultichain::new(),
            }
        );

        assert_eq!(op.dependencies(), &deps);

        assert_eq!(
            op1.hash_key(),
            format!(
                "BasisTimeEvolution({}; [{}]; {}; [{}])",
                false,
                deps.hash_key(),
                braket.hash_key(),
                PertMultichain::new().hash_key(),
            )
        );
        assert!(!op1.is_scalar());
        assert_eq!(format!("{}", op1), "op(T)");

        let op2 = BasisTimeEvolution::builder(make_super_multichain(&deps, 1u32)).build().unwrap();

        assert_ne!(&op1, &op2);
    }

    #[test]
    fn test_differentiation() {
        let len_pert_name: u32 = 2;
        let deps = make_pert_multichain(len_pert_name, 8u32, 1u32, 10u32);
        let op = BasisTimeEvolution::builder(deps.clone()).build().unwrap();

        let p: Arc<Perturbation> = deps.keys().first().cloned().unwrap();
        let mut diff_op = op.differentiate(&p).unwrap();
        let mut deriv = PertMultichain::new();
        deriv.insert(&p);

        let diff_cast = downcast_from_arc::<BasisTimeEvolution>(&diff_op).unwrap();

        assert_eq!(diff_cast.derivative(), &deriv);

        let max_order = deps.get_order(&p);
        for _ in 1..=2 * max_order + 1 {
            diff_op = diff_op.differentiate(&p).unwrap();
        }

        assert_eq!(&diff_op, &ZeroOperator::new());

        assert!(is_zero_expr(
            &op.differentiate(&make_perturbation_symbol(len_pert_name + 1u32, 4u32)).unwrap(),
            None
        ));
    }

    #[test]
    fn test_serialization() {
        let op = BasisTimeEvolution::builder(make_pert_multichain(2u32, 8u32, 1u32, 10u32))
            .build()
            .unwrap();
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let deps = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
        let op1 = BasisTimeEvolution::builder(deps.clone()).build().unwrap();

        assert!(is_expr_type::<BasisTimeEvolution>(&op1));
        assert!(!is_zero_expr(&op1, None));
        assert!(!is_one_expr(&op1, None));

        let op2 = BasisTimeEvolution::builder(deps.clone()).build().unwrap();
        let op3 = BasisTimeEvolution::builder(make_super_multichain(&deps, 1u32)).build().unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));
    }
}
