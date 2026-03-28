use std::collections::{BTreeMap, HashSet};
use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::internal::{any_interned_expr_matches, intern_expr, join_mapped};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{
    NumberTolerance, downcast_from_arc, expression_error, generic_expression_error,
    get_number_tolerance, is_zero_expr, unreachable_error,
};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct EliminationRule {
    parameter: Arc<dyn Expr>,
    min_order: u32,
    perturbations: Vec<Arc<Perturbation>>,
}

impl EliminationRule {
    pub fn new(
        parameter: Arc<dyn Expr>,
        min_order: u32,
        perturbations: &[Arc<Perturbation>],
    ) -> Self {
        let mut perturbations = perturbations.to_vec();
        perturbations.sort();

        Self {
            parameter,
            min_order,
            perturbations,
        }
    }

    #[inline]
    pub fn parameter(&self) -> &Arc<dyn Expr> {
        &self.parameter
    }

    #[inline]
    pub fn min_order(&self) -> u32 {
        self.min_order
    }

    #[inline]
    pub fn perturbations(&self) -> &[Arc<Perturbation>] {
        &self.perturbations
    }

    #[inline]
    pub fn hash_key(&self) -> String {
        format!(
            "EliminationRule({}; {}; [{}])",
            self.parameter.hash_key(),
            self.min_order,
            join_mapped(self.perturbations.iter(), ",", |p| p.hash_key()),
        )
    }
}

impl PartialEq for EliminationRule {
    fn eq(&self, other: &Self) -> bool {
        &self.parameter == &other.parameter
            && self.min_order == other.min_order
            && self.perturbations == other.perturbations
    }
}

impl Eq for EliminationRule {}

impl std::fmt::Display for EliminationRule {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "({}; {}; [{}])",
            self.parameter,
            self.min_order,
            join_mapped(self.perturbations.iter(), ",", |p| p.to_string()),
        )
    }
}

// A `SubExpr` is uniquely determined by its fields `name`, `expression` and `derivative`.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct SubExpr {
    name: String,
    expression: Arc<dyn Expr>,
    derivative: PertMultichain,
    // The use of elimination rules is mostly to let users track which
    // parameters have been eliminated. We do not use it for comparison.
    elimination_rules: Vec<EliminationRule>,
    // `zero_rules_applied` is mostly used by the function `apply_zero_rules()`. We
    // do not use it for equality comparison, either.
    zero_rules_applied: bool,
}

impl SubExpr {
    #[inline]
    pub fn builder(name: impl Into<String>, expression: Arc<dyn Expr>) -> SubExprBuilder {
        SubExprBuilder {
            name: name.into(),
            expression,
            derivative: None,
            elimination_rules: None,
            zero_rules_applied: None,
            // We will always check the name conflict when users try to build a `SubExpr`
            check_name_conflict: true,
        }
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[inline]
    pub fn expression(&self) -> &Arc<dyn Expr> {
        &self.expression
    }

    #[inline]
    pub fn derivative(&self) -> &PertMultichain {
        &self.derivative
    }

    #[inline]
    pub fn elimination_rules(&self) -> &[EliminationRule] {
        &self.elimination_rules
    }

    #[inline]
    pub fn zero_rules_applied(&self) -> bool {
        self.zero_rules_applied
    }
}

// This function is used for checking whether users try to build a new
// undifferentiated `SubExpr` with the same name but different `expression`
#[inline]
fn subexpr_name_conflict(name: &str, expression: &Arc<dyn Expr>) -> bool {
    any_interned_expr_matches(|expr| {
        let Some(subexpr) = downcast_from_arc::<SubExpr>(expr) else {
            return false;
        };

        subexpr.name == name
            && subexpr.derivative.is_empty()
            && subexpr.elimination_rules.is_empty()
            && &subexpr.expression != expression
    })
}

#[derive(Debug)]
pub struct SubExprBuilder {
    name: String,
    expression: Arc<dyn Expr>,
    derivative: Option<PertMultichain>,
    elimination_rules: Option<Vec<EliminationRule>>,
    zero_rules_applied: Option<bool>,
    check_name_conflict: bool,
}

impl SubExprBuilder {
    // The following methods are ONLY used inside this file
    #[inline]
    fn derivative(mut self, derivative: PertMultichain) -> Self {
        self.derivative = Some(derivative);
        self
    }

    #[inline]
    fn elimination_rules(mut self, elimination_rules: Vec<EliminationRule>) -> Self {
        self.elimination_rules = Some(elimination_rules);
        self
    }

    #[inline]
    fn zero_rules_applied(mut self, zero_rules_applied: bool) -> Self {
        self.zero_rules_applied = Some(zero_rules_applied);
        self
    }

    // This function should MOSTLY be called by `SubExpr::builder()` function
    #[inline]
    fn check_name_conflict(mut self, check_name_conflict: bool) -> Self {
        self.check_name_conflict = check_name_conflict;
        self
    }

    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        if is_zero_expr(&self.expression, None) {
            return impl_zero_expr!(self.expression.is_scalar());
        }

        let derivative = self.derivative.unwrap_or(PertMultichain::new());
        if self.check_name_conflict && !derivative.is_empty() {
            return Err(unreachable_error(
                format!(
                    "A new SubExpr {} cannot have derivative {}",
                    self.name,
                    derivative.to_string()
                ),
                &self.expression,
                None,
            ));
        }

        let mut elimination_rules = self.elimination_rules.unwrap_or_default();
        elimination_rules.sort_by_key(|rule| rule.hash_key());
        if self.check_name_conflict && !elimination_rules.is_empty() {
            return Err(unreachable_error(
                format!(
                    "A new SubExpr {} cannot have elimination rules {}",
                    self.name,
                    join_mapped(elimination_rules.iter(), ";", |rule| rule.to_string()),
                ),
                &self.expression,
                None,
            ));
        }

        let zero_rules_applied = self.zero_rules_applied.unwrap_or(false);

        // Check whether users try to build a new undifferentiated `SubExpr`
        // with the same name but different `expression`
        if self.check_name_conflict {
            if subexpr_name_conflict(&self.name, &self.expression) {
                return Err(expression_error(
                    format!(
                        "A SubExpr exists with the same name {} but different expression",
                        self.name
                    ),
                    &self.expression,
                    None,
                ));
            }
        }

        Ok(intern_expr(Arc::new(SubExpr {
            name: self.name,
            expression: self.expression,
            derivative,
            elimination_rules,
            zero_rules_applied,
        })))
    }
}

impl ExprInternal for SubExpr {
    // It is more appropriate to set `zero_rules_applied` as false after the
    // functions `replace_expr_children`, `retain_expr_fields` and `replace_expr_self`
    impl_unary_expr_internal_methods!(
        SubExpr,
        Argument,
        expression,
        true,
        |this: &SubExpr, arg| {
            SubExpr::builder(this.name.clone(), arg)
                .derivative(this.derivative.clone())
                .elimination_rules(this.elimination_rules.clone())
                .check_name_conflict(false)
                .build()
        }
    );

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "SubExpr({}; {{{}}}; [{}]; [{}]; {})",
            self.name,
            self.expression.hash_key(),
            self.derivative.hash_key(),
            join_mapped(self.elimination_rules.iter(), ";", |rule| rule.hash_key()),
            self.zero_rules_applied,
        )
    }

    #[inline]
    fn total_order(&self) -> u32 {
        self.derivative.total_order()
    }

    // For functions `deep_eq_superchains()` and `eq_by_superchains()`, our
    // policy is that we ignore the defail of `SubExpr`, i.e. we ignore the
    // field `expression`, but treat `SubExpr` as a free symbol.
    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<SubExpr>(other) {
            // Since this function is used by `find_superchains()`, we check ONLY `name` and `derivative`.
            self.name == op.name && self.derivative.is_subchain(&op.derivative)
        } else {
            false
        }
    }

    #[inline]
    fn eq_by_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<SubExpr>(other) {
            // This function is used by `replace()` and `retain_expr()`, except
            // for `name` and `derivative`, we also require the same value of
            // `zero_rules_applied`.
            self.name == op.name
                && self.derivative.is_subchain(&op.derivative)
                && self.zero_rules_applied == op.zero_rules_applied
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for SubExpr {
    #[inline]
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        self.expression.is_scalar()
    }

    #[inline]
    fn clone_expr(&self) -> Arc<dyn Expr> {
        Arc::new(self.clone())
    }

    #[inline]
    fn apply_zero_rules(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.zero_rules_applied {
            return Ok(self.clone_expr());
        }

        let new_expr = self.expression.apply_zero_rules(freq_tol.clone()).map_err(|e| {
            generic_expression_error(
                format!(
                    "SubExpr::apply_zero_rules() failed with tolerance {}",
                    freq_tol.clone().unwrap_or_else(get_number_tolerance)
                ),
                self,
                Some(Box::new(e)),
            )
        })?;

        if is_zero_expr(&new_expr, freq_tol) {
            return impl_zero_expr!(new_expr.is_scalar());
        } else {
            SubExpr::builder(self.name.clone(), new_expr)
                .derivative(self.derivative.clone())
                .elimination_rules(self.elimination_rules.clone())
                .zero_rules_applied(true)
                .check_name_conflict(false)
                .build()
        }
    }

    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_expr = self.expression.differentiate(s).map_err(|e| {
            generic_expression_error(
                format!(
                    "SubExpr::differentiate() failed for differentiation with respect to {}",
                    s
                ),
                self,
                Some(Box::new(e)),
            )
        })?;

        SubExpr::builder(self.name.clone(), diff_expr)
            .derivative(self.derivative.with_added_perturbation(s))
            .elimination_rules(self.elimination_rules.clone())
            .zero_rules_applied(self.zero_rules_applied)
            .check_name_conflict(false)
            .build()
    }

    #[inline]
    fn eliminate(
        &self,
        parameter: &Arc<dyn Expr>,
        perturbations: &[Arc<Perturbation>],
        min_order: u32,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.elimination_rules.iter().any(|rule| rule.parameter() == parameter) {
            return Err(generic_expression_error(
                format!(
                    "SubExpr::eliminate() got repeated elimination of a parameter {}",
                    parameter
                ),
                self,
                None,
            ));
        }

        let new_expr =
            self.expression.eliminate(parameter, perturbations, min_order).map_err(|e| {
                generic_expression_error(
                    format!(
                        "SubExpr::eliminate() failed for the parameter {}, minimum order {}, and perturbations [{}]",
                        parameter,
                        min_order,
                        join_mapped(perturbations, ",", |p| p.to_string())
                    ),
                    self,
                    Some(Box::new(e)),
                )
            })?;

        let mut elimination_rules = self.elimination_rules.clone();
        elimination_rules.push(EliminationRule::new(parameter.clone(), min_order, perturbations));

        SubExpr::builder(self.name.clone(), new_expr)
            .derivative(self.derivative.clone())
            .elimination_rules(elimination_rules)
            .zero_rules_applied(self.zero_rules_applied)
            .check_name_conflict(false)
            .build()
    }

    #[inline]
    fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>, include_derivatives: bool) -> bool {
        self.match_self_any(set, include_derivatives)
            || self.expression.exist_any(set, include_derivatives)
    }

    #[inline]
    fn find_superchains(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
        if self.deep_eq_superchains(s) {
            BTreeMap::from([(self.total_order(), HashSet::from([self.clone_expr()]))])
        } else {
            self.expression.find_superchains(s)
        }
    }

    #[inline]
    fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_self_any(set, false) {
            return impl_zero_expr!(self.expression.is_scalar());
        }

        let new_expr = self.expression.remove(set).map_err(|e| {
            generic_expression_error(
                format!(
                    "SubExpr::remove() failed for removing {{{}}}",
                    join_mapped(set.iter(), ",", |s| s.to_string()),
                ),
                self,
                Some(Box::new(e)),
            )
        })?;

        if &new_expr == &self.expression {
            Ok(self.clone_expr())
        } else {
            SubExpr::builder(self.name.clone(), new_expr)
                .derivative(self.derivative.clone())
                .elimination_rules(self.elimination_rules.clone())
                .zero_rules_applied(self.zero_rules_applied)
                .check_name_conflict(false)
                .build()
        }
    }

    #[inline]
    fn retain(
        &self,
        set: &HashSet<Arc<dyn Expr>>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_self_any(set, include_derivatives) {
            return Ok(self.clone_expr());
        }

        let new_expr = self.expression.retain(set, include_derivatives).map_err(|e| {
            generic_expression_error(
                format!(
                    "SubExpr::retain() failed for set {{{}}}",
                    join_mapped(set.iter(), ",", |s| s.to_string()),
                ),
                self,
                Some(Box::new(e)),
            )
        })?;

        if &new_expr == &self.expression {
            Ok(self.clone_expr())
        } else {
            SubExpr::builder(self.name.clone(), new_expr)
                .derivative(self.derivative.clone())
                .elimination_rules(self.elimination_rules.clone())
                .zero_rules_applied(self.zero_rules_applied)
                .check_name_conflict(false)
                .build()
        }
    }
}

//FIXME: should we allow the comparison between SubExpr and other Expr in the function `eq_expr()`?
impl PartialEq for SubExpr {
    fn eq(&self, other: &Self) -> bool {
        // We also compare `expression`, which may change after some methods
        // like `apply_zero_rules()`, `remove()`, `replace()` and `retain()`.
        self.name == other.name
            && &self.expression == &other.expression
            && self.derivative == other.derivative
        //&& self.elimination_rules == other.elimination_rules
        //&& self.zero_rules_applied == other.zero_rules_applied
    }
}

impl Eq for SubExpr {}

impl std::fmt::Display for SubExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.elimination_rules.is_empty() {
            write!(
                f,
                "{}({}; {{{}}})^{}",
                self.name, self.zero_rules_applied, self.expression, self.derivative
            )
        } else {
            write!(
                f,
                "{}([{}]; {}; {{{}}})^{}",
                self.name,
                join_mapped(self.elimination_rules.iter(), ",", |rule| rule.to_string()),
                self.zero_rules_applied,
                self.expression,
                self.derivative,
            )
        }
    }
}
