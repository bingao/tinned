use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::expressions::{OneElecOperator, WfnParameter, ZeroOperator};
use crate::perturbations::Perturbation;
use crate::utils::{intern, invalid_expression_error, is_zero_expr};

/// A TemporumOperator is a non-scalar operator acting on a ket or a bra
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct TemporumOperator {
    on_ket: bool,
    argument: Arc<dyn Expr>,
}

impl TemporumOperator {
    #[inline]
    pub fn builder(argument: Arc<dyn Expr>) -> TemporumOperatorBuilder {
        TemporumOperatorBuilder { on_ket: true, argument }
    }

    #[inline]
    fn builder_from(&self, argument: Arc<dyn Expr>) -> TemporumOperatorBuilder {
        TemporumOperatorBuilder { on_ket: self.on_ket, argument }
    }

    #[inline]
    pub fn on_ket(&self) -> bool {
        self.on_ket
    }

    #[inline]
    pub fn argument(&self) -> &Arc<dyn Expr> {
        &self.argument
    }
}

#[derive(Debug)]
pub struct TemporumOperatorBuilder {
    on_ket: bool,
    argument: Arc<dyn Expr>,
}

impl TemporumOperatorBuilder {
    #[inline]
    pub fn on_ket(mut self, on_ket: bool) -> Self {
        self.on_ket = on_ket;
        self
    }

    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.argument.is_scalar() {
            return Err(invalid_expression_error(
                "TemporumOperatorBuilder::build() - scalar argument",
                &self.argument,
            ));
        }

        if self.argument.is::<OneElecOperator>() || self.argument.is::<WfnParameter>() {
            Ok(intern(Arc::new(TemporumOperator { on_ket: self.on_ket, argument: self.argument })))
        } else {
            Err(invalid_expression_error(
                "TemporumOperatorBuilder::build() - unsupported argument type",
                &self.argument,
            ))
        }
    }
}

impl Expr for TemporumOperator {
    #[inline]
    fn as_any(&self) -> &dyn Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        format!("TemporumOperator({}; {})", self.on_ket, self.argument.hash_key())
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        false
    }

    fn differentiate(&self, s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_arg = self.argument.differentiate(s)?;

        if is_zero_expr(&diff_arg) {
            Ok(ZeroOperator::new())
        } else {
            self.builder_from(diff_arg).build()
        }
    }
}

impl Display for TemporumOperator {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "{}({})", if self.on_ket { "i*dt" } else { "-i*dt" }, self.argument)
    }
}
