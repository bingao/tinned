use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::internal::{intern_expr, join_mapped};
use crate::perturbations::{PertMultichain, Perturbation};
use crate::public::{
    NumberTolerance, downcast_from_arc, generic_expression_error, get_number_tolerance,
    is_zero_expr,
};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct SubExpr {
    name: String,
    expression: Arc<dyn Expr>,
    derivative: PertMultichain,
    elimination_rules: HashMap<Arc<dyn Expr>, (u32, Vec<Arc<Perturbation>>)>,
    is_zero_strength: bool,
}

impl SubExpr {
    #[inline]
    pub fn builder(name: impl Into<String>, expression: Arc<dyn Expr>) -> SubExprBuilder {
        SubExprBuilder {
            name: name.into(),
            expression,
            derivative: None,
            elimination_rules: None,
            is_zero_strength: None,
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
    pub fn elimination_rules(&self) -> &HashMap<Arc<dyn Expr>, (u32, Vec<Arc<Perturbation>>)> {
        &self.elimination_rules
    }

    #[inline]
    pub fn is_zero_strength(&self) -> bool {
        self.is_zero_strength
    }
}

#[derive(Debug)]
pub struct SubExprBuilder {
    name: String,
    expression: Arc<dyn Expr>,
    derivative: Option<PertMultichain>,
    elimination_rules: Option<HashMap<Arc<dyn Expr>, (u32, Vec<Arc<Perturbation>>)>>,
    is_zero_strength: Option<bool>,
}

impl SubExprBuilder {
    #[inline]
    fn derivative(mut self, derivative: PertMultichain) -> Self {
        self.derivative = Some(derivative);
        self
    }

    #[inline]
    fn elimination_rules(
        mut self,
        elimination_rules: HashMap<Arc<dyn Expr>, (u32, Vec<Arc<Perturbation>>)>,
    ) -> Self {
        self.elimination_rules = Some(elimination_rules);
        self
    }

    #[inline]
    fn is_zero_strength(mut self, is_zero_strength: bool) -> Self {
        self.is_zero_strength = Some(is_zero_strength);
        self
    }

    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        if is_zero_expr(&self.expression, None) {
            return impl_zero_expr!(self.expression.is_scalar());
        }

        let derivative = self.derivative.unwrap_or(PertMultichain::new());
        let elimination_rules = self.elimination_rules.unwrap_or(HashMap::new());
        let is_zero_strength = self.is_zero_strength.unwrap_or(false);

        Ok(intern_expr(Arc::new(SubExpr {
            name: self.name,
            expression: self.expression,
            derivative,
            elimination_rules,
            is_zero_strength,
        })))
    }
}

impl ExprInternal for SubExpr {
    // It is more appropriate to set `is_zero_strength` as false after the
    // functions `replace_expr_fields`, `retain_expr_fields` and `replace_expr_self`
    impl_unary_expr_internal_methods!(SubExpr, expression, true, |this: &SubExpr, arg| {
        SubExpr::builder(this.name.clone(), arg)
            .derivative(this.derivative.clone())
            .elimination_rules(this.elimination_rules.clone())
            .build()
    });

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "SubExpr({}; {}; {}; [{}]; {})",
            self.name,
            self.expression.hash_key(),
            self.derivative.hash_key(),
            join_mapped(
                self.elimination_rules.iter(),
                ",",
                |(parameter, (min_order, perturbations))| {
                    format!(
                        "{}([{}]; {})",
                        parameter.hash_key(),
                        join_mapped(perturbations, ",", |p| p.hash_key()),
                        min_order
                    )
                }
            ),
            self.is_zero_strength,
        )
    }

    #[inline]
    fn total_order(&self) -> u32 {
        self.derivative.total_order()
    }

    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<SubExpr>(other) {
            // We find sub expressions with
            // `is_zero_strength` either `true` or `false`
            self.name == op.name
                && self.expression.deep_eq_superchains(&op.expression)
                && self.derivative.is_subchain(&op.derivative)
                && self.elimination_rules == op.elimination_rules
        } else {
            false
        }
    }

    #[inline]
    fn eq_by_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<SubExpr>(other) {
            self.name == op.name
                && self.expression.deep_eq_superchains(&op.expression)
                && self.derivative.is_subchain(&op.derivative)
                && self.elimination_rules == op.elimination_rules
                && self.is_zero_strength == op.is_zero_strength
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
    fn clean_temporum(
        &self,
        freq_tol: Option<NumberTolerance>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.is_zero_strength {
            return Ok(self.clone_expr());
        }

        let new_expr = self.expression.clean_temporum(freq_tol.clone()).map_err(|e| {
            generic_expression_error(
                format!(
                    "SubExpr::clean_temporum() failed with tolerance {}",
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
                .is_zero_strength(true)
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
            .is_zero_strength(self.is_zero_strength)
            .build()
    }

    #[inline]
    fn eliminate(
        &self,
        parameter: &Arc<dyn Expr>,
        perturbations: &[Arc<Perturbation>],
        min_order: u32,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.elimination_rules.contains_key(parameter) {
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
        elimination_rules.insert(parameter.clone(), (min_order, perturbations.to_vec()));

        SubExpr::builder(self.name.clone(), new_expr)
            .derivative(self.derivative.clone())
            .elimination_rules(elimination_rules)
            .is_zero_strength(self.is_zero_strength)
            .build()
    }

    #[inline]
    fn exist_any(&self, set: &HashSet<Arc<dyn Expr>>) -> bool {
        set.iter().any(|expr| self.eq_expr(expr.as_ref())) || self.expression.exist_any(set)
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
        if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
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
                .is_zero_strength(self.is_zero_strength)
                .build()
        }
    }
}

impl PartialEq for SubExpr {
    fn eq(&self, other: &Self) -> bool {
        // We also compare `expression`, which may change after `clean_temporum()`
        self.name == other.name
            && &self.expression == &other.expression
            && self.derivative == other.derivative
            && self.elimination_rules == other.elimination_rules
            && self.is_zero_strength == other.is_zero_strength
    }
}

impl Eq for SubExpr {}

impl std::fmt::Display for SubExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.elimination_rules.is_empty() {
            write!(f, "{}({})^{}", self.name, self.is_zero_strength, self.derivative)
        } else {
            write!(
                f,
                "{}([{}]; {})^{}",
                self.name,
                join_mapped(
                    self.elimination_rules.iter(),
                    ",",
                    |(parameter, (min_order, perturbations))| {
                        format!(
                            "{}([{}]; {})",
                            parameter,
                            join_mapped(perturbations, ",", |p| p.to_string()),
                            min_order
                        )
                    }
                ),
                self.is_zero_strength,
                self.derivative,
            )
        }
    }
}
