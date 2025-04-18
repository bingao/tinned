use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{HermitianTranspose, Number, Transpose, ZeroOperator};
use crate::utils::is_expr_type;

/// Dot product of a bra and a ket (inner product)
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct DotProduct {
    bra: Arc<dyn Expr>,
    ket: Arc<dyn Expr>,
    is_complex: bool,
}

impl DotProduct {
    pub fn new(
        bra: Arc<dyn Expr>,
        ket: Arc<dyn Expr>,
        is_complex: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if bra.is_scalar() || ket.is_scalar() {
            return Err(crate::utils::invalid_expression_error(
                "DotProduct::new() - both arguments must be non-scalar",
                if bra.is_scalar() {
                    &bra
                } else {
                    &ket
                },
            ));
        }

        if is_expr_type::<ZeroOperator>(&bra) || is_expr_type::<ZeroOperator>(&ket) {
            return Ok(Number::zero());
        }

        let bra = if is_complex {
            HermitianTranspose::new(bra)?
        } else {
            Transpose::new(bra)?
        };

        Ok(crate::utils::intern_expr(Arc::new(Self {
            bra,
            ket,
            is_complex,
        })))
    }

    #[inline]
    pub fn bra(&self) -> &Arc<dyn Expr> {
        &self.bra
    }

    #[inline]
    pub fn ket(&self) -> &Arc<dyn Expr> {
        &self.ket
    }

    #[inline]
    pub fn is_complex(&self) -> bool {
        self.is_complex
    }
}

#[typetag::serde]
impl Expr for DotProduct {
    #[inline]
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        format!("DotProduct({}; {}; {})", self.bra.hash_key(), self.ket.hash_key(), self.is_complex)
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        true
    }

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        if let Some(dot) = crate::utils::downcast_from_ref::<DotProduct>(other) {
            self == dot
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
        let diff_bra = self.bra.differentiate(s)?;
        let diff_ket = self.ket.differentiate(s)?;

        Self::new(diff_bra, diff_ket, self.is_complex)
    }
}

impl PartialEq for DotProduct {
    fn eq(&self, other: &Self) -> bool {
        self.is_complex == other.is_complex && &self.bra == &other.bra && &self.ket == &other.ket
    }
}

impl Eq for DotProduct {}

impl std::fmt::Display for DotProduct {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "<{}, {}>", self.bra, self.ket)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_struct_safety!(DotProduct);

    test_thread_interning!({
        DotProduct::new(make_wfn_parameter("bra"), make_wfn_parameter("ket"), true).unwrap()
    });

    #[test]
    fn test_impl_expr() {
        let psi1 = make_wfn_parameter("");
        let mut op0 = DotProduct::new(ZeroOperator::new(), psi1.clone(), true).unwrap();
        assert!(is_zero_expr(&op0));

        op0 = DotProduct::new(psi1.clone(), ZeroOperator::new(), false).unwrap();
        assert!(is_zero_expr(&op0));

        let psi1_dagger = HermitianTranspose::new(psi1.clone()).unwrap();
        let psi2 = make_wfn_parameter("");
        let op1 = DotProduct::new(psi1.clone(), psi2.clone(), true).unwrap();

        assert_eq!(
            op1.hash_key(),
            format!("DotProduct({}; {}; {})", psi1_dagger.hash_key(), psi2.hash_key(), true)
        );
        assert!(op1.is_scalar());
        assert_eq!(format!("{}", op1), format!("<{}, {}>", psi1_dagger, psi2));

        let mut op = downcast_from_arc::<DotProduct>(&op1).unwrap();
        assert_eq!(
            op,
            &DotProduct {
                bra: psi1_dagger.clone(),
                ket: psi2.clone(),
                is_complex: true
            }
        );

        assert_eq!(op.bra(), &psi1_dagger);
        assert_eq!(op.ket(), &psi2);
        assert!(op.is_complex());

        let op2 = DotProduct::new(psi1.clone(), psi2.clone(), true).unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        let op3 = DotProduct::new(psi1.clone(), psi2.clone(), false).unwrap();
        assert_ne!(&op1, &op3);

        op = downcast_from_arc::<DotProduct>(&op3).unwrap();
        let psi1_tr = Transpose::new(psi1.clone()).unwrap();
        assert_eq!(
            op,
            &DotProduct {
                bra: psi1_tr.clone(),
                ket: psi2.clone(),
                is_complex: false
            }
        );
        assert_eq!(op.bra(), &psi1_tr);
        assert_eq!(op.ket(), &psi2);
        assert!(!op.is_complex());

        let op4 = DotProduct::new(psi2, psi1.clone(), true).unwrap();
        assert_ne!(&op1, &op4);

        let op5 = DotProduct::new(psi1.clone(), psi1.clone(), true).unwrap();
        assert_ne!(&op1, &op5);
    }

    #[test]
    fn test_serialization() {
        let op = DotProduct::new(make_wfn_parameter(""), make_wfn_parameter(""), true).unwrap();
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let psi1 = make_wfn_parameter("");
        let psi2 = make_wfn_parameter("");
        let op1 = DotProduct::new(psi1.clone(), psi2.clone(), true).unwrap();
        let op2 = DotProduct::new(psi1.clone(), psi2.clone(), true).unwrap();
        let op3 = DotProduct::new(psi1.clone(), psi2.clone(), false).unwrap();
        let op4 = DotProduct::new(psi2.clone(), psi1.clone(), true).unwrap();
        let op5 = DotProduct::new(psi2.clone(), psi1.clone(), false).unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));
        assert!(!Arc::ptr_eq(&op1, &op4));
        assert!(!Arc::ptr_eq(&op1, &op5));

        assert!(is_expr_type::<DotProduct>(&op1));
        assert!(!is_zero_expr(&op1));
        assert!(!is_one_expr(&op1));
    }
}
