use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::expressions::{
    Add, HermitianTranspose, MatrixMul, Mul, Number, Power, Transpose, ZeroOperator,
};
use crate::perturbations::Perturbation;
use crate::utils::{downcast_expr, intern, is_one_expr};

/// Represents complex conjugation of an expression.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Conjugate {
    argument: Arc<dyn Expr>,
}

impl Conjugate {
    pub fn new(argument: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if let Some(num) = downcast_expr::<Number>(&argument) {
            return Ok(num.conjugate().into());
        } else if let Some(add) = downcast_expr::<Add>(&argument) {
            let terms: Vec<Arc<dyn Expr>> =
                add.terms().iter().map(|t| Self::new(t.clone())?).collect();
            return Add::new(terms);
        } else if let Some(mul) = downcast_expr::<Mul>(&argument) {
            let coef = mul.coefficient().conjugate();
            let mut new_terms: Vec<Arc<dyn Expr>> =
                mul.factors().iter().map(|f| Self::new(f.clone())?).collect();
            new_terms.push(coef.into());
            return Mul::new(new_terms);
        } else if let Some(power) = downcast_expr::<Power>(&argument) {
            let new_base = Self::new(power.base().clone())?;
            return Power::new(new_base, power.exponent());
        } else if argument.is::<ZeroOperator>() {
            return Ok(argument);
        } else if let Some(conj) = downcast_expr::<Conjugate>(&argument) {
            return Ok(conj.argument.clone());
        } else if let Some(trans) = downcast_expr::<Transpose>(&argument) {
            return HermitianTranspose::new(trans.argument().clone());
        } else if let Some(herm) = downcast_expr::<HermitianTranspose>(&argument) {
            return Transpose::new(herm.argument().clone());
        } else if let Some(matmul) = downcast_expr::<MatrixMul>(&argument) {
            if is_one_expr(matmul.coefficient()) {
                return Ok(intern(Arc::new(Self { argument })));
            }

            let new_arg = MatrixMul::new(matmul.factors().to_vec())?;
            return MatrixMul::new(vec![
                Self::new(matmul.coefficient().clone())?,
                intern(Arc::new(Self { argument: new_arg })),
            ]);
        }

        Ok(intern(Arc::new(Self { argument })))
    }

    #[inline]
    pub fn argument(&self) -> &Arc<dyn Expr> {
        &self.argument
    }
}

impl Expr for Conjugate {
    #[inline]
    fn as_any(&self) -> &dyn Any { self }

    #[inline]
    fn hash_key(&self) -> String {
        format!("Conjugate({})", self.argument.hash_key())
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        self.argument.is_scalar()
    }

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        if let Some(conj) = downcast_expr::<Conjugate>(other) {
            self.argument == conj.argument
        } else {
            false
        }
    }

    fn differentiate(&self, s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_arg = self.argument.differentiate(s)?;
        Self::new(diff_arg)
    }
}

impl Display for Conjugate {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "conj({arg})", arg = self.argument)
    }
}
