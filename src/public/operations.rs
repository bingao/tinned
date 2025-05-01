use std::sync::Arc;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Add, MatrixAdd, MatrixMul, Mul, Number, Power};
use crate::perturbations::Perturbation;
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
pub fn differentiate_expr(
    expr: &Arc<dyn Expr>,
    perturbations: &[Arc<Perturbation>],
) -> Result<Arc<dyn Expr>, TinnedError> {
    let mut deriv = expr.clone();
    for pert in perturbations {
        deriv = deriv
            .differentiate(pert)
            .map_err(|e| expression_error("Differentiation failed", expr, Some(Box::new(e))))?;
    }

    Ok(deriv)
}
