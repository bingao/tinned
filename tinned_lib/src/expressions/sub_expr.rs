use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;
use std::sync::Arc;

use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::internal::{intern_expr, join_mapped};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{
    NumberTolerance, downcast_from_arc, generic_expression_error, get_number_tolerance,
    is_zero_expr,
};

#[derive(Clone, Debug, Serialize, Deserialize)]
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
}

impl fmt::Display for EliminationRule {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "({}; {}; [{}])",
            self.parameter,
            self.min_order,
            join_mapped(self.perturbations.iter(), ",", |p| p.to_string()),
        )
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ReplacementRule {
    #[serde(with = "replacement_rule_map_serde")]
    map: HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
    include_derivatives: bool,
}

impl ReplacementRule {
    pub fn new(map: HashMap<Arc<dyn Expr>, Arc<dyn Expr>>, include_derivatives: bool) -> Self {
        Self {
            map,
            include_derivatives,
        }
    }

    #[inline]
    pub fn map(&self) -> &HashMap<Arc<dyn Expr>, Arc<dyn Expr>> {
        &self.map
    }

    #[inline]
    pub fn include_derivatives(&self) -> bool {
        self.include_derivatives
    }
}

impl fmt::Display for ReplacementRule {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{{")?;

        for (k, v) in &self.map {
            write!(f, "{} -> {}; ", k, v)?;
        }

        write!(f, "{}}}", self.include_derivatives)
    }
}

mod replacement_rule_map_serde {
    use super::*;

    pub fn serialize<S>(
        map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        serializer: S,
    ) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let entries: Vec<(&Arc<dyn Expr>, &Arc<dyn Expr>)> = map.iter().collect();

        entries.serialize(serializer)
    }

    pub fn deserialize<'de, D>(
        deserializer: D,
    ) -> Result<HashMap<Arc<dyn Expr>, Arc<dyn Expr>>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let entries: Vec<(Arc<dyn Expr>, Arc<dyn Expr>)> = Vec::deserialize(deserializer)?;

        Ok(entries.into_iter().collect())
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RetainmentRule {
    s: Arc<dyn Expr>,
    include_derivatives: bool,
}

impl RetainmentRule {
    pub fn new(s: Arc<dyn Expr>, include_derivatives: bool) -> Self {
        Self {
            s,
            include_derivatives,
        }
    }

    #[inline]
    pub fn s(&self) -> &Arc<dyn Expr> {
        &self.s
    }

    #[inline]
    pub fn include_derivatives(&self) -> bool {
        self.include_derivatives
    }
}

impl fmt::Display for RetainmentRule {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{{{}; {}}}", self.s, self.include_derivatives)
    }
}

// A `SubExpr` represents one high level concrete expression struct with
// `name`, and `expression` containing its detail.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SubExpr {
    name: String,
    expression: Arc<dyn Expr>,
    // Users may give one same `name` but different `expression`'s when they
    // build `SubExpr`, which is not reasonable. We do not prevent it, but we
    // introduce a unique `identifier` that combines `name` and the hash key of
    // `expression` for the first time the `SubExpr` was built. Instead of
    // using `name`, we actually use `identifier` to distinguish one `SubExpr`
    // from others. Note that `identifier` should not change in different
    // methods of the trait `Expr`.
    identifier: String,
    derivative: PertMultichain,
    // The following fields are mostly used to let users track which methods
    // have been applied to `SubExpr`. We do not use them for comparison of
    // `SubExpr`'s.
    at_zero_perturbations: bool,
    elimination_rules: Vec<EliminationRule>,
    removal_rules: Vec<HashSet<Arc<dyn Expr>>>,
    replacement_rules: Vec<ReplacementRule>,
    retainment_rules: Vec<RetainmentRule>,
}

impl SubExpr {
    #[inline]
    pub fn new(name: impl Into<String>, expression: Arc<dyn Expr>) -> Arc<dyn Expr> {
        let name = name.into();
        let identifier = format!("{}({})", name, expression.hash_value());

        intern_expr(Arc::new(Self {
            name,
            expression,
            identifier,
            derivative: PertMultichain::new(),
            at_zero_perturbations: false,
            elimination_rules: Vec::new(),
            removal_rules: Vec::new(),
            replacement_rules: Vec::new(),
            retainment_rules: Vec::new(),
        }))
    }

    #[inline]
    fn with_zero_perturbations(&self, expression: Arc<dyn Expr>) -> Arc<dyn Expr> {
        intern_expr(Arc::new(Self {
            name: self.name.clone(),
            expression,
            identifier: self.identifier.clone(),
            derivative: self.derivative.clone(),
            at_zero_perturbations: true,
            elimination_rules: self.elimination_rules.clone(),
            removal_rules: self.removal_rules.clone(),
            replacement_rules: self.replacement_rules.clone(),
            retainment_rules: self.retainment_rules.clone(),
        }))
    }

    #[inline]
    fn with_differentiation(
        &self,
        expression: Arc<dyn Expr>,
        s: Arc<Perturbation>,
    ) -> Arc<dyn Expr> {
        intern_expr(Arc::new(Self {
            name: self.name.clone(),
            expression,
            identifier: self.identifier.clone(),
            derivative: self.derivative.with_added_perturbation(s),
            at_zero_perturbations: self.at_zero_perturbations,
            elimination_rules: self.elimination_rules.clone(),
            removal_rules: self.removal_rules.clone(),
            replacement_rules: self.replacement_rules.clone(),
            retainment_rules: self.retainment_rules.clone(),
        }))
    }

    #[inline]
    fn with_elimination(
        &self,
        expression: Arc<dyn Expr>,
        parameter: Arc<dyn Expr>,
        perturbations: &[Arc<Perturbation>],
        min_order: u32,
    ) -> Arc<dyn Expr> {
        let mut elimination_rules = self.elimination_rules.clone();
        elimination_rules.push(EliminationRule::new(parameter, min_order, perturbations));

        intern_expr(Arc::new(Self {
            name: self.name.clone(),
            expression,
            identifier: self.identifier.clone(),
            derivative: self.derivative.clone(),
            at_zero_perturbations: self.at_zero_perturbations,
            elimination_rules,
            removal_rules: self.removal_rules.clone(),
            replacement_rules: self.replacement_rules.clone(),
            retainment_rules: self.retainment_rules.clone(),
        }))
    }

    #[inline]
    fn with_removal(
        &self,
        expression: Arc<dyn Expr>,
        set: HashSet<Arc<dyn Expr>>,
    ) -> Arc<dyn Expr> {
        let mut removal_rules = self.removal_rules.clone();
        removal_rules.push(set);

        intern_expr(Arc::new(Self {
            name: self.name.clone(),
            expression,
            identifier: self.identifier.clone(),
            derivative: self.derivative.clone(),
            at_zero_perturbations: self.at_zero_perturbations,
            elimination_rules: self.elimination_rules.clone(),
            removal_rules,
            replacement_rules: self.replacement_rules.clone(),
            retainment_rules: self.retainment_rules.clone(),
        }))
    }

    #[inline]
    fn with_replacement(
        &self,
        expression: Arc<dyn Expr>,
        map: HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        include_derivatives: bool,
    ) -> Arc<dyn Expr> {
        let mut replacement_rules = self.replacement_rules.clone();
        replacement_rules.push(ReplacementRule::new(map, include_derivatives));

        intern_expr(Arc::new(Self {
            name: self.name.clone(),
            expression,
            identifier: self.identifier.clone(),
            derivative: self.derivative.clone(),
            at_zero_perturbations: self.at_zero_perturbations,
            elimination_rules: self.elimination_rules.clone(),
            removal_rules: self.removal_rules.clone(),
            replacement_rules,
            retainment_rules: self.retainment_rules.clone(),
        }))
    }

    #[inline]
    fn with_retainment(
        &self,
        expression: Arc<dyn Expr>,
        s: Arc<dyn Expr>,
        include_derivatives: bool,
    ) -> Arc<dyn Expr> {
        let mut retainment_rules = self.retainment_rules.clone();
        retainment_rules.push(RetainmentRule::new(s, include_derivatives));

        intern_expr(Arc::new(Self {
            name: self.name.clone(),
            expression,
            identifier: self.identifier.clone(),
            derivative: self.derivative.clone(),
            at_zero_perturbations: self.at_zero_perturbations,
            elimination_rules: self.elimination_rules.clone(),
            removal_rules: self.removal_rules.clone(),
            replacement_rules: self.replacement_rules.clone(),
            retainment_rules,
        }))
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
    pub fn identifier(&self) -> &str {
        &self.identifier
    }

    #[inline]
    pub fn derivative(&self) -> &PertMultichain {
        &self.derivative
    }

    #[inline]
    pub fn at_zero_perturbations(&self) -> bool {
        self.at_zero_perturbations
    }

    #[inline]
    pub fn elimination_rules(&self) -> &[EliminationRule] {
        &self.elimination_rules
    }

    #[inline]
    pub fn removal_rules(&self) -> &[HashSet<Arc<dyn Expr>>] {
        &self.removal_rules
    }

    #[inline]
    pub fn replacement_rules(&self) -> &[ReplacementRule] {
        &self.replacement_rules
    }

    #[inline]
    pub fn retainment_rules(&self) -> &[RetainmentRule] {
        &self.retainment_rules
    }
}

impl ExprInternal for SubExpr {
    impl_expr_internal_methods!(SubExpr, true);

    #[inline]
    fn replace_one_in_children(
        &self,
        expr: &Arc<dyn Expr>,
        replacement: Arc<dyn Expr>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        impl_unary_expr_arg_operation!(
            self,
            Argument,
            expression,
            |arg: &Arc<dyn Expr>| arg.replace_one(expr, replacement.clone(), include_derivatives),
            "SubExpr::replace_one_in_children() failed",
            |this: &SubExpr, arg| Ok(this.with_replacement(
                arg,
                HashMap::from([(expr.clone(), replacement)]),
                include_derivatives
            ))
        )
    }

    #[inline]
    fn replace_all_in_children(
        &self,
        map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        impl_unary_expr_arg_operation!(
            self,
            Argument,
            expression,
            |arg: &Arc<dyn Expr>| arg.replace_all(map, include_derivatives),
            "SubExpr::replace_all_in_children() failed",
            |this: &SubExpr, arg| Ok(this.with_replacement(arg, map.clone(), include_derivatives))
        )
    }

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "SubExpr({}; {{{}}}; [{}])",
            self.identifier,
            self.expression.hash_key(),
            self.derivative.hash_key(),
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
            // Since this function is used by `find_all()`, we check
            // ONLY `identifier` and `derivative`.
            self.identifier == op.identifier && self.derivative.is_subchain(&op.derivative)
        } else {
            false
        }
    }

    #[inline]
    fn eq_by_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<SubExpr>(other) {
            self.identifier == op.identifier && self.derivative.is_subchain(&op.derivative)
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
    fn clone_expr(&self) -> Arc<dyn Expr> {
        intern_expr(Arc::new(self.clone()))
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        self.expression.is_scalar()
    }

    #[inline]
    fn has_unperturbed_term(&self) -> bool {
        self.expression.has_unperturbed_term()
    }

    #[inline]
    fn substitute_zero_perturbations(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let new_expr =
            self.expression.substitute_zero_perturbations(freq_tol.clone()).map_err(|e| {
                generic_expression_error(
                    format!(
                        "SubExpr::substitute_zero_perturbations() failed with tolerance {}",
                        freq_tol.clone().unwrap_or_else(get_number_tolerance)
                    ),
                    self,
                    Some(Box::new(e)),
                )
            })?;

        if is_zero_expr(&new_expr, freq_tol) {
            impl_zero_expr!(new_expr.is_scalar())
        } else {
            Ok(self.with_zero_perturbations(new_expr))
        }
    }

    fn differentiate(&self, s: Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_expr = self.expression.differentiate(s.clone()).map_err(|e| {
            generic_expression_error(
                format!(
                    "SubExpr::differentiate() failed for differentiation with respect to {}",
                    s
                ),
                self,
                Some(Box::new(e)),
            )
        })?;

        Ok(self.with_differentiation(diff_expr, s))
    }

    #[inline]
    fn eliminate(
        &self,
        parameter: Arc<dyn Expr>,
        perturbations: &[Arc<Perturbation>],
        min_order: u32,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.elimination_rules.iter().any(|rule| rule.parameter() == &parameter) {
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
            self.expression.eliminate(parameter.clone(), perturbations, min_order).map_err(|e| {
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

        Ok(self.with_elimination(new_expr, parameter, perturbations, min_order))
    }

    #[inline]
    fn find_all(&self, s: &Arc<dyn Expr>) -> BTreeMap<u32, HashSet<Arc<dyn Expr>>> {
        if self.deep_eq_superchains(s) {
            BTreeMap::from([(self.total_order(), HashSet::from([self.clone_expr()]))])
        } else {
            self.expression.find_all(s)
        }
    }

    #[inline]
    fn match_one(&self, s: &Arc<dyn Expr>, include_derivatives: bool) -> bool {
        self.match_one_self(s, include_derivatives)
            || self.expression.match_one(s, include_derivatives)
    }

    #[inline]
    fn match_any(&self, set: &HashSet<Arc<dyn Expr>>, include_derivatives: bool) -> bool {
        self.match_any_self(set, include_derivatives)
            || self.expression.match_any(set, include_derivatives)
    }

    #[inline]
    fn remove_one(&self, s: &Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_one_self(s, false) {
            return impl_zero_expr!(self.expression.is_scalar());
        }

        let new_expr = self.expression.remove_one(s).map_err(|e| {
            generic_expression_error(
                format!("SubExpr::remove_one() failed for removing {}", s,),
                self,
                Some(Box::new(e)),
            )
        })?;

        if &new_expr == &self.expression {
            Ok(self.clone_expr())
        } else {
            Ok(self.with_removal(new_expr, HashSet::from([s.clone()])))
        }
    }

    #[inline]
    fn remove_all(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_any_self(set, false) {
            return impl_zero_expr!(self.expression.is_scalar());
        }

        let new_expr = self.expression.remove_all(set).map_err(|e| {
            generic_expression_error(
                format!(
                    "SubExpr::remove_all() failed for removing {{{}}}",
                    join_mapped(set.iter(), ",", |s| s.to_string()),
                ),
                self,
                Some(Box::new(e)),
            )
        })?;

        if &new_expr == &self.expression {
            Ok(self.clone_expr())
        } else {
            Ok(self.with_removal(new_expr, set.clone()))
        }
    }

    #[inline]
    fn retain_one(
        &self,
        s: &Arc<dyn Expr>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_one_self(s, include_derivatives) {
            return Ok(self.clone_expr());
        }

        impl_unary_expr_arg_operation!(
            self,
            Argument,
            expression,
            |arg: &Arc<dyn Expr>| arg.retain_one(s, include_derivatives),
            "SubExpr::retain_one() failed",
            |this: &SubExpr, arg| Ok(this.with_retainment(arg, s.clone(), include_derivatives))
        )
    }
}

//FIXME: should we allow the comparison between SubExpr and other Expr in the function `eq_expr()`?
impl PartialEq for SubExpr {
    fn eq(&self, other: &Self) -> bool {
        // We also compare `expression`, which may change after some methods
        // like `substitute_zero_perturbations()`, `remove_all()`, `replace_all()` and `retain_one()`.
        self.identifier == other.identifier
            && &self.expression == &other.expression
            && self.derivative == other.derivative
    }
}

impl Eq for SubExpr {}

impl fmt::Display for SubExpr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}({{{}}}; {}; [{}]; [{}]; [{}]; [{}])^{}",
            self.identifier,
            self.expression,
            self.at_zero_perturbations,
            join_mapped(self.elimination_rules.iter(), ";", |rule| rule.to_string()),
            join_mapped(self.removal_rules.iter(), ";", |exprs| join_mapped(
                exprs.iter(),
                ";",
                |expr| expr.to_string()
            )),
            join_mapped(self.replacement_rules.iter(), ";", |rule| rule.to_string()),
            join_mapped(self.retainment_rules.iter(), ";", |rule| rule.to_string()),
            self.derivative,
        )
    }
}
