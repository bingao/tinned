use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{
    Add, HermitianTranspose, MatrixMul, Mul, Number, Power, Transpose, ZeroOperator,
};
use crate::utils::{downcast_from_arc, downcast_from_ref, intern_expr, is_expr_type, is_one_expr};

/// Represents complex conjugation of an expression.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Conjugate {
    argument: Arc<dyn Expr>,
}

impl Conjugate {
    pub fn new(argument: Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if let Some(num) = downcast_from_arc::<Number>(&argument) {
            return Ok(num.conjugate().into());
        } else if let Some(add) = downcast_from_arc::<Add>(&argument) {
            let terms: Vec<Arc<dyn Expr>> =
                add.terms().iter().map(|t| Self::new(t.clone())).collect::<Result<_, _>>()?;
            return Add::new(terms);
        } else if let Some(mul) = downcast_from_arc::<Mul>(&argument) {
            let coef = mul.coefficient().conjugate();
            let mut new_terms: Vec<Arc<dyn Expr>> =
                mul.factors().iter().map(|f| Self::new(f.clone())).collect::<Result<_, _>>()?;
            new_terms.push(coef.into());
            return Mul::new(new_terms);
        } else if let Some(power) = downcast_from_arc::<Power>(&argument) {
            let new_base = Self::new(power.base().clone())?;
            return Power::new(new_base, power.exponent());
        } else if is_expr_type::<ZeroOperator>(&argument) {
            return Ok(argument);
        } else if let Some(conj) = downcast_from_arc::<Conjugate>(&argument) {
            return Ok(conj.argument.clone());
        } else if let Some(trans) = downcast_from_arc::<Transpose>(&argument) {
            return HermitianTranspose::new(trans.argument().clone());
        } else if let Some(herm) = downcast_from_arc::<HermitianTranspose>(&argument) {
            return Transpose::new(herm.argument().clone());
        } else if let Some(matmul) = downcast_from_arc::<MatrixMul>(&argument) {
            if is_one_expr(matmul.coefficient()) {
                return Ok(intern_expr(Arc::new(Self { argument })));
            }

            let new_arg = MatrixMul::new(matmul.factors().to_vec())?;
            return MatrixMul::new(vec![
                Self::new(matmul.coefficient().clone())?,
                intern_expr(Arc::new(Self { argument: new_arg })),
            ]);
        }

        Ok(intern_expr(Arc::new(Self { argument })))
    }

    #[inline]
    pub fn argument(&self) -> &Arc<dyn Expr> {
        &self.argument
    }
}

#[typetag::serde]
impl Expr for Conjugate {
    #[inline]
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

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
        if let Some(conj) = downcast_from_ref::<Conjugate>(other) {
            self == conj
        } else {
            false
        }
    }

    #[inline]
    fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }

    fn differentiate(
        &self,
        s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_arg = self.argument.differentiate(s)?;
        Self::new(diff_arg)
    }
}

impl PartialEq for Conjugate {
    fn eq(&self, other: &Self) -> bool {
        &self.argument == &other.argument
    }
}

impl Eq for Conjugate {}

impl std::fmt::Display for Conjugate {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "conj({arg})", arg = self.argument)
    }
}
