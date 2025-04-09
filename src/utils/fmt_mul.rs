use std::fmt::{Formatter, Result as FmtResult};

use crate::expressions::{MatrixMul, Mul};
use crate::utils::is_one_expr;

// Helper function of stream formatting of Mul
// Format: coef * factor1 * factor2 * ..., omit coef if one
#[inline]
pub fn fmt_mul(f: &mut Formatter, mul: &Mul) -> FmtResult {
    let mut wrote_any = false;

    let coef = mul.coefficient();
    if coef.is_one() {
        write!(f, "{}", coef)?;
        wrote_any = true;
    }

    for factor in mul.factors() {
        if wrote_any {
            write!(f, " * ")?;
        }
        write!(f, "{}", factor)?;
        wrote_any = true;
    }

    Ok(())
}

// Helper function of stream formatting of MatrixMul
// Format: coef * factor1 * factor2 * ..., omit coef if one
#[inline]
pub fn fmt_matrix_mul(f: &mut Formatter, mul: &MatrixMul) -> FmtResult {
    let mut wrote_any = false;

    let coef = mul.coefficient();
    if !is_one_expr(coef) {
        write!(f, "{}", coef)?;
        wrote_any = true;
    }

    for factor in mul.factors() {
        if wrote_any {
            write!(f, " * ")?;
        }
        write!(f, "{}", factor)?;
        wrote_any = true;
    }

    Ok(())
}
