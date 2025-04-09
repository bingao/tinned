// Helper function of stream formatting of Mul
// Format: coef * factor1 * factor2 * ..., omit coef if one
#[inline]
pub fn fmt_mul(f: &mut std::fmt::Formatter, mul: &crate::expressions::Mul) -> std::fmt::Result {
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
pub fn fmt_matrix_mul(
    f: &mut std::fmt::Formatter,
    mul: &crate::expressions::MatrixMul,
) -> std::fmt::Result {
    let mut wrote_any = false;

    let coef = mul.coefficient();
    if !crate::utils::is_one_expr(coef) {
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
