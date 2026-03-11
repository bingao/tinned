use std::error::Error;
use std::sync::Arc;

use crate::core::{Expr, TinnedError};
use crate::internal::join_mapped;
use crate::perturbations::Perturbation;

#[inline]
pub fn multi_expression_error(
    message: impl Into<String>,
    exprs: &[Arc<dyn Expr>],
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    TinnedError::ExpressionError {
        message: message.into(),
        expression: join_mapped(exprs, ";", |expr| expr.to_string()),
        source,
    }
}

#[inline]
pub fn multi_perturbation_error(
    message: impl Into<String>,
    perturbations: &[Arc<Perturbation>],
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    TinnedError::PerturbationError {
        message: message.into(),
        perturbation: join_mapped(perturbations, ";", |p| p.to_string()),
        source,
    }
}

#[inline]
pub fn expression_error(
    message: impl Into<String>,
    expr: &Arc<dyn Expr>,
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    multi_expression_error(message, std::slice::from_ref(expr), source)
}

#[inline]
pub fn generic_expression_error<E: Expr>(
    message: impl Into<String>,
    expr: &E,
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    TinnedError::ExpressionError {
        message: message.into(),
        expression: format!("{:?}", expr),
        source,
    }
}

#[inline]
pub fn perturbation_error(
    message: impl Into<String>,
    perturbation: &Arc<Perturbation>,
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    multi_perturbation_error(message, std::slice::from_ref(perturbation), source)
}

#[inline]
pub fn unreachable_error(
    message: impl Into<String>,
    expr: &Arc<dyn Expr>,
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    TinnedError::Unreachable {
        message: message.into(),
        expression: format!("{:?}", expr),
        source,
    }
}

#[inline]
pub fn generic_error(
    message: impl Into<String>,
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    TinnedError::GenericError {
        message: message.into(),
        source,
    }
}
