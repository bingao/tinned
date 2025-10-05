use std::error::Error;
use std::sync::Arc;

use crate::core::{Expr, TinnedError};
use crate::internal::{multi_expression_format, multi_perturbation_format};
use crate::perturbations::Perturbation;

#[inline]
pub fn multi_expression_error(
    message: &'static str,
    exprs: &[Arc<dyn Expr>],
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    TinnedError::ExpressionError {
        message,
        expression: multi_expression_format(exprs, ";"),
        source,
    }
}

#[inline]
pub fn multi_perturbation_error(
    message: &'static str,
    perts: &[Arc<Perturbation>],
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    TinnedError::PerturbationError {
        message,
        perturbation: multi_perturbation_format(perts, ";"),
        source,
    }
}

#[inline]
pub fn expression_error(
    message: &'static str,
    expr: &Arc<dyn Expr>,
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    multi_expression_error(message, std::slice::from_ref(expr), source)
}

#[inline]
pub fn generic_expression_error<E: Expr>(
    message: &'static str,
    expr: &E,
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    TinnedError::ExpressionError {
        message,
        expression: format!("{:?}", expr),
        source,
    }
}

#[inline]
pub fn perturbation_error(
    message: &'static str,
    pert: &Arc<Perturbation>,
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    multi_perturbation_error(message, std::slice::from_ref(pert), source)
}

#[inline]
pub fn unreachable_error(
    message: &'static str,
    expr: &Arc<dyn Expr>,
    source: Option<Box<dyn Error + Send + Sync>>,
) -> TinnedError {
    TinnedError::Unreachable {
        message,
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
