use std::sync::Arc;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Add, MatrixAdd, MatrixMul, Mul, Number, Power};
use crate::perturbations::PertSequence;
use crate::public::{downcast_from_arc, expression_error, multi_expression_error};

#[inline]
pub fn negate_expr(expr: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
    if let Some(num) = downcast_from_arc::<Number>(&expr) {
        return Ok(num.negate().into());
    }

    let minus_one = Number::minus_one();

    if expr.is_scalar() {
        Mul::new(vec![minus_one, expr])
    } else {
        MatrixMul::new(vec![minus_one, expr])
    }
}

#[inline]
pub fn subtract_exprs(
    lhs: Arc<dyn Expr>,
    rhs: Arc<dyn Expr>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    match (lhs.is_scalar(), rhs.is_scalar()) {
        (true, true) => Add::new(vec![lhs, negate_expr(rhs)?]),
        (false, false) => MatrixAdd::new(vec![lhs, negate_expr(rhs)?]),
        _ => Err(multi_expression_error(
            "subtract_exprs: lhs and rhs must both be scalar or both be non-scalar",
            &vec![lhs, rhs],
            None,
        )),
    }
}

#[inline]
pub fn divide_exprs(
    numerator: Arc<dyn Expr>,
    denominator: Arc<dyn Expr>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    Mul::new(vec![numerator, Power::new(denominator, -1)?])
}

/// Performs high-order differentiation
#[inline]
pub fn differentiate_expr<T: PertSequence>(
    expr: &Arc<dyn Expr>,
    perturbations: &T,
) -> Result<Arc<dyn Expr>, TinnedError> {
    let mut result = expr.clone();
    for pert in perturbations.ordered_perturbations() {
        result = result
            .differentiate(&pert)
            .map_err(|e| expression_error("Differentiation failed", expr, Some(Box::new(e))))?;
    }

    Ok(result)
}

/// Commutator [a, b] = ab - ba
#[inline]
pub fn commutator(a: Arc<dyn Expr>, b: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
    if a.is_scalar() {
        return Err(expression_error("Commutator requires non-scalar a", &a, None));
    }
    if b.is_scalar() {
        return Err(expression_error("Commutator requires non-scalar b", &b, None));
    }

    subtract_exprs(MatrixMul::new(vec![a.clone(), b.clone()])?, MatrixMul::new(vec![b, a])?)
}

/// Anticommutator {a, b} = ab + ba
#[inline]
pub fn anticommutator(a: Arc<dyn Expr>, b: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
    if a.is_scalar() {
        return Err(expression_error("Anticommutator requires non-scalar a", &a, None));
    }
    if b.is_scalar() {
        return Err(expression_error("Anticommutator requires non-scalar b", &b, None));
    }

    MatrixAdd::new(vec![MatrixMul::new(vec![a.clone(), b.clone()])?, MatrixMul::new(vec![b, a])?])
}

/// S commutator [a, b]_{s} = asb - bsa
#[inline]
pub fn s_commutator(
    a: Arc<dyn Expr>,
    b: Arc<dyn Expr>,
    s: Arc<dyn Expr>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    if a.is_scalar() {
        return Err(expression_error("S commutator requires non-scalar a", &a, None));
    }
    if b.is_scalar() {
        return Err(expression_error("S commutator requires non-scalar b", &b, None));
    }
    if s.is_scalar() {
        return Err(expression_error("S commutator requires non-scalar s", &s, None));
    }

    subtract_exprs(
        MatrixMul::new(vec![a.clone(), s.clone(), b.clone()])?,
        MatrixMul::new(vec![b, s, a])?,
    )
}

/// S anticommutator {a, b}_{s} = asb + bsa
#[inline]
pub fn s_anticommutator(
    a: Arc<dyn Expr>,
    b: Arc<dyn Expr>,
    s: Arc<dyn Expr>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    if a.is_scalar() {
        return Err(expression_error("S anticommutator requires non-scalar a", &a, None));
    }
    if b.is_scalar() {
        return Err(expression_error("S anticommutator requires non-scalar b", &b, None));
    }
    if s.is_scalar() {
        return Err(expression_error("S anticommutator requires non-scalar s", &s, None));
    }

    MatrixAdd::new(vec![
        MatrixMul::new(vec![a.clone(), s.clone(), b.clone()])?,
        MatrixMul::new(vec![b, s, a])?,
    ])
}

/// Returns sum of perturbations' frequencies
#[inline]
pub fn sum_pert_frequencies<T: PertSequence>(
    perturbations: &T,
) -> Result<Arc<dyn Expr>, TinnedError> {
    Add::new(
        perturbations.ordered_perturbations().into_iter().map(|p| p.frequency().clone()).collect(),
    )
}
