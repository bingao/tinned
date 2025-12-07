use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{
    Add, Conjugate, HermitianTranspose, MatrixMul, Mul, Number, Transpose, ZeroOperator,
};
use crate::internal::intern_expr;
use crate::perturbations::Perturbation;
use crate::public::{
    NumberTolerance, downcast_from_arc, downcast_from_ref, expression_error,
    generic_expression_error, is_expr_type, is_one_expr,
};

/// Dot product of a bra and a ket (inner product)
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct DotProduct {
    bra: Arc<dyn Expr>,
    ket: Arc<dyn Expr>,
    allow_braket_swap: bool,
    is_scalar: bool,
}

impl DotProduct {
    pub fn new(
        bra: Arc<dyn Expr>,
        use_hermitian: bool,
        ket: Arc<dyn Expr>,
        allow_braket_swap: bool,
        is_scalar: Option<bool>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let bra = if use_hermitian {
            HermitianTranspose::new(bra)?
        } else {
            Transpose::new(bra)?
        };

        let is_scalar = is_scalar.unwrap_or(true);

        Self::make_dot_product(bra, ket, allow_braket_swap, is_scalar)
    }

    // Helper function to strip scalar coefficient from a MatrixMul expression
    #[inline]
    fn strip_matrixmul_coefficient(
        expr: Arc<dyn Expr>,
        coefficients: &mut Vec<Arc<dyn Expr>>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if let Some(mat_mul) = downcast_from_arc::<MatrixMul>(&expr) {
            if !is_one_expr(mat_mul.coefficient(), None) {
                coefficients.push(mat_mul.coefficient().clone());
                return MatrixMul::new(mat_mul.factors().to_vec());
            }
        }
        Ok(expr)
    }

    #[inline]
    fn make_dot_product(
        bra: Arc<dyn Expr>,
        ket: Arc<dyn Expr>,
        allow_braket_swap: bool,
        is_scalar: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if bra.is_scalar() {
            return Err(expression_error("DotProduct requires bra must be non-scalar", &bra, None));
        }
        if ket.is_scalar() {
            return Err(expression_error("DotProduct requires ket must be non-scalar", &ket, None));
        }

        if is_expr_type::<ZeroOperator>(&bra) || is_expr_type::<ZeroOperator>(&ket) {
            return Ok(Number::zero());
        }

        let mut coefficients = Vec::with_capacity(2);

        let bra = Self::strip_matrixmul_coefficient(bra, &mut coefficients)?;
        let ket = Self::strip_matrixmul_coefficient(ket, &mut coefficients)?;

        let (bra, ket) = if allow_braket_swap {
            let trans_ket = Transpose::new(ket.clone())?;
            if bra.hash_key() > trans_ket.hash_key() {
                (trans_ket, Transpose::new(bra)?)
            } else {
                (bra, ket)
            }
        } else {
            (bra, ket)
        };

        let product = intern_expr(Arc::new(Self {
            bra,
            ket,
            allow_braket_swap,
            is_scalar,
        }));

        if coefficients.is_empty() {
            Ok(product)
        } else {
            coefficients.push(product);
            Mul::new(coefficients)
        }
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
    pub fn allow_braket_swap(&self) -> bool {
        self.allow_braket_swap
    }

    #[inline]
    pub fn conjugate(&self) -> Result<Arc<dyn Expr>, TinnedError> {
        impl_binary_expr_arg_operation!(
            self,
            bra,
            ket,
            |arg: &Arc<dyn Expr>| Conjugate::new(arg.clone()),
            "DotProduct::conjugate() failed",
            |this: &DotProduct, bra, ket| Self::make_dot_product(bra, ket, this.allow_braket_swap, this.is_scalar)
        )
    }
}

impl ExprInternal for DotProduct {
    impl_binary_expr_internal_methods!(
        DotProduct,
        bra,
        ket,
        false,
        |this: &DotProduct, bra, ket| Self::make_dot_product(bra, ket, this.allow_braket_swap, this.is_scalar)
    );

    #[inline]
    fn hash_key(&self) -> String {
        // Important to include `allow_braket_swap` into hash
        format!(
            "DotProduct({}; {}; {})",
            self.bra.hash_key(),
            self.ket.hash_key(),
            self.allow_braket_swap,
        )
    }

    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(dot) = downcast_from_arc::<DotProduct>(other) {
            if self.bra.deep_eq_superchains(&dot.bra) && self.ket.deep_eq_superchains(&dot.ket) {
                return true;
            }

            if self.allow_braket_swap || dot.allow_braket_swap {
                let trans_bra = match Transpose::new(self.bra.clone()) {
                    Ok(expr) => expr,
                    Err(_) => return false,
                };
                if !trans_bra.deep_eq_superchains(&dot.ket) {
                    return false;
                }
                let trans_ket = match Transpose::new(self.ket.clone()) {
                    Ok(expr) => expr,
                    Err(_) => return false,
                };

                return trans_ket.deep_eq_superchains(&dot.bra);
            }

            false
        } else {
            false
        }
    }

    // For unambiguous replacement, we require equality of `bra` and `ket`, and
    // make replacement by ignoring derivatives of `bra` and `ket`. So we do
    // not override methods `eq_by_superchains()` and `replace_expr_self()` of
    // `ExprInternal`.
}

#[typetag::serde]
impl Expr for DotProduct {
    impl_binary_expr_common_methods!(
        DotProduct,
        bra,
        ket,
        is_scalar,
        |this: &DotProduct, bra, ket| Self::make_dot_product(bra, ket, this.allow_braket_swap, this.is_scalar),
        true
    );

    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_bra = self.bra.differentiate(s).map_err(|e| {
            generic_expression_error(
                "DotProduct::differentiate() failed for bra",
                self,
                Some(Box::new(e)),
            )
        })?;
        let diff_ket = self.ket.differentiate(s).map_err(|e| {
            generic_expression_error(
                "DotProduct::differentiate() failed for ket",
                self,
                Some(Box::new(e)),
            )
        })?;

        Add::new(vec![
            Self::make_dot_product(diff_bra, self.ket.clone(), self.allow_braket_swap, true)?,
            Self::make_dot_product(self.bra.clone(), diff_ket, self.allow_braket_swap, true)?,
        ])
    }
}

impl PartialEq for DotProduct {
    fn eq(&self, other: &Self) -> bool {
        if &self.bra == &other.bra && &self.ket == &other.ket {
            return true;
        }

        if self.allow_braket_swap || other.allow_braket_swap {
            let trans_bra = match Transpose::new(self.bra.clone()) {
                Ok(expr) => expr,
                Err(_) => return false,
            };
            if &trans_bra != &other.ket {
                return false;
            }
            let trans_ket = match Transpose::new(self.ket.clone()) {
                Ok(expr) => expr,
                Err(_) => return false,
            };

            return &trans_ket == &other.bra;
        }

        false
    }
}

impl Eq for DotProduct {}

impl std::fmt::Display for DotProduct {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.allow_braket_swap {
            write!(f, "<{};{}>", self.bra, self.ket)
        } else {
            write!(f, "<{}|{}>", self.bra, self.ket)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::symbol::test_utils::make_symbol;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::is_zero_expr;

    test_struct_safety!(DotProduct);

    test_thread_interning!({
        DotProduct::new(make_wfn_parameter("bra"), true, make_wfn_parameter("ket"), true, None).unwrap()
    });

    #[test]
    fn test_impl_expr() {
        let allow_braket_swap = true;
        let use_hermitian = true;
        let psi1 = make_wfn_parameter("");
        let mut op0 =
            DotProduct::new(ZeroOperator::new(), use_hermitian, psi1.clone(), allow_braket_swap, None)
                .unwrap();
        assert!(is_zero_expr(&op0, None));

        op0 = DotProduct::new(psi1.clone(), !use_hermitian, ZeroOperator::new(), allow_braket_swap, None)
            .unwrap();
        assert!(is_zero_expr(&op0, None));

        let psi2 = make_wfn_parameter("");
        let op1 =
            DotProduct::new(psi1.clone(), use_hermitian, psi2.clone(), allow_braket_swap, None).unwrap();

        let psi1_dagger = HermitianTranspose::new(psi1.clone()).unwrap();
        let conj_psi1 = Conjugate::new(psi1.clone()).unwrap();
        let trans_psi2 = Transpose::new(psi2.clone()).unwrap();

        let (bra, ket) = if psi1_dagger.hash_key() <= trans_psi2.hash_key() {
            (psi1_dagger, psi2.clone())
        } else {
            (trans_psi2.clone(), conj_psi1)
        };

        assert_eq!(
            op1.hash_key(),
            format!("DotProduct({}; {}; {})", bra.hash_key(), ket.hash_key(), allow_braket_swap)
        );
        assert!(op1.is_scalar());
        assert_eq!(format!("{}", op1), format!("<{};{}>", bra, ket));

        let mut op = downcast_from_arc::<DotProduct>(&op1).unwrap();
        assert_eq!(
            op,
            &DotProduct {
                bra: bra.clone(),
                ket: ket.clone(),
                allow_braket_swap,
                is_scalar: true,
            }
        );

        assert_eq!(op.bra(), &bra);
        assert_eq!(op.ket(), &ket);
        assert!(op.allow_braket_swap());
        assert_eq!(
            &op.conjugate().unwrap(),
            &DotProduct::new(
                ket.clone(),
                true,
                HermitianTranspose::new(bra.clone()).unwrap(),
                true,
                None,
            )
            .unwrap()
        );

        let op2 =
            DotProduct::new(psi1.clone(), use_hermitian, psi2.clone(), allow_braket_swap, None).unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        let op3 =
            DotProduct::new(psi1.clone(), !use_hermitian, psi2.clone(), allow_braket_swap, None).unwrap();
        assert_ne!(&op1, &op3);

        let op4 =
            DotProduct::new(psi2.clone(), use_hermitian, psi1.clone(), allow_braket_swap, None).unwrap();
        assert_ne!(&op1, &op4);

        let op5 =
            DotProduct::new(psi1.clone(), use_hermitian, psi1.clone(), allow_braket_swap, None).unwrap();
        assert_ne!(&op1, &op5);

        let op6 =
            DotProduct::new(psi2.clone(), !use_hermitian, psi1.clone(), allow_braket_swap, None).unwrap();
        assert!(Arc::ptr_eq(&op3, &op6));
        assert_eq!(&op3, &op6);

        let op7 = DotProduct::new(psi2.clone(), !use_hermitian, psi1.clone(), !allow_braket_swap, None)
            .unwrap();
        assert!(!Arc::ptr_eq(&op3, &op7));
        assert_eq!(&op3, &op7);
        assert_eq!(format!("{}", op7), format!("<{}|{}>", trans_psi2, psi1));

        op = downcast_from_arc::<DotProduct>(&op7).unwrap();
        assert_eq!(
            op,
            &DotProduct {
                bra: trans_psi2.clone(),
                ket: psi1.clone(),
                allow_braket_swap: !allow_braket_swap,
                is_scalar: true
            }
        );

        assert_eq!(op.bra(), &trans_psi2);
        assert_eq!(op.ket(), &psi1);
        assert!(!op.allow_braket_swap());
        assert_eq!(
            &op.conjugate().unwrap(),
            &DotProduct::new(psi2.clone(), true, Conjugate::new(psi1.clone()).unwrap(), false, None)
                .unwrap()
        );

        let coef_psi1 = make_symbol(4u32);
        let coef_psi2 = make_symbol(4u32);
        let op8 = DotProduct::new(
            MatrixMul::new(vec![coef_psi1.clone(), psi1.clone()]).unwrap(),
            use_hermitian,
            MatrixMul::new(vec![coef_psi2.clone(), psi2.clone()]).unwrap(),
            allow_braket_swap,
            None,
        )
        .unwrap();

        assert!(is_expr_type::<Mul>(&op8));
        assert_eq!(
            &op8,
            &Mul::new(vec![
                Conjugate::new(coef_psi1.clone()).unwrap(),
                coef_psi2.clone(),
                op1.clone(),
            ])
            .unwrap()
        );
    }

    #[test]
    fn test_differentiation() {
        let psi1 = make_wfn_parameter("");
        let psi2 = make_wfn_parameter("");
        let op = DotProduct::new(psi1.clone(), true, psi2.clone(), true, None).unwrap();
        let p = make_perturbation_symbol(4u32, 4u32);
        let mut diff_op = op.differentiate(&p).unwrap();
        let diff_psi1 = psi1.differentiate(&p).unwrap();
        let diff_psi2 = psi2.differentiate(&p).unwrap();

        assert_eq!(
            &diff_op,
            &Add::new(vec![
                DotProduct::new(diff_psi1.clone(), true, psi2.clone(), true, None).unwrap(),
                DotProduct::new(psi1.clone(), true, diff_psi2.clone(), true, None).unwrap(),
            ])
            .unwrap()
        );

        diff_op = diff_op.differentiate(&p).unwrap();

        assert_eq!(
            &diff_op,
            &Add::new(vec![
                DotProduct::new(diff_psi1.differentiate(&p).unwrap(), true, psi2.clone(), true, None)
                    .unwrap(),
                Mul::new(vec![
                    Number::from_i64(2),
                    DotProduct::new(diff_psi1.clone(), true, diff_psi2.clone(), true, None).unwrap()
                ])
                .unwrap(),
                DotProduct::new(psi1.clone(), true, diff_psi2.differentiate(&p).unwrap(), true, None)
                    .unwrap(),
            ])
            .unwrap()
        );
    }

    #[test]
    fn test_serialization() {
        let op =
            DotProduct::new(make_wfn_parameter(""), true, make_wfn_parameter(""), true, None).unwrap();
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let psi1 = make_wfn_parameter("");
        let psi2 = make_wfn_parameter("");
        let op = DotProduct::new(psi1.clone(), true, psi2.clone(), true, None).unwrap();

        assert!(is_expr_type::<DotProduct>(&op));
        assert!(!is_zero_expr(&op, None));
        assert!(!is_one_expr(&op, None));

        let op1 = DotProduct::new(psi1.clone(), true, psi2.clone(), true, None).unwrap();
        let op2 = DotProduct::new(psi1.clone(), true, psi2.clone(), false, None).unwrap();
        let op3 = DotProduct::new(psi1.clone(), false, psi2.clone(), true, None).unwrap();
        let op4 = DotProduct::new(psi1.clone(), false, psi2.clone(), false, None).unwrap();

        assert!(Arc::ptr_eq(&op, &op1));
        assert!(!Arc::ptr_eq(&op, &op2));
        assert!(!Arc::ptr_eq(&op, &op3));
        assert!(!Arc::ptr_eq(&op, &op4));

        assert_eq!(&op, &op2);
        assert_eq!(&op3, &op4);

        let op5 = DotProduct::new(psi2.clone(), true, psi1.clone(), true, None).unwrap();
        let op6 = DotProduct::new(psi2.clone(), true, psi1.clone(), false, None).unwrap();
        let op7 = DotProduct::new(psi2.clone(), false, psi1.clone(), true, None).unwrap();
        let op8 = DotProduct::new(psi2.clone(), false, psi1.clone(), false, None).unwrap();

        assert!(!Arc::ptr_eq(&op, &op5));
        assert!(!Arc::ptr_eq(&op, &op6));
        assert!(!Arc::ptr_eq(&op, &op7));
        assert!(!Arc::ptr_eq(&op, &op8));

        assert_eq!(&op5, &op6);
        assert_eq!(&op3, &op7);
        assert_eq!(&op3, &op8);
    }
}
