use std::fmt::{Formatter, Result as FmtResult};
use std::sync::Arc;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Composition, MatrixMul, Mul, Trace, WfnParameter};
use crate::utils::{invalid_expression_error, is_one_expr};

// Helper function to validate the density matrix, grid weight and overlap
// distribution
#[inline]
pub fn validate_xc_inputs(
    density_matrix: &Arc<dyn Expr>,
    grid_weight: &Arc<dyn Expr>,
    overlap_distribution: &Arc<dyn Expr>,
) -> Result<(), TinnedError> {
    if !density_matrix.is::<WfnParameter>() {
        return Err(invalid_expression_error(
            "validate_xc_inputs() - density matrix must be WfnParameter",
            density_matrix,
        ));
    }

    if !grid_weight.is_scalar() {
        return Err(invalid_expression_error(
            "validate_xc_inputs() - grid weight must be scalar",
            grid_weight,
        ));
    }

    if overlap_distribution.is_scalar() {
        return Err(invalid_expression_error(
            "validate_xc_inputs() - overlap distribution must be non-scalar",
            overlap_distribution,
        ));
    }

    Ok(())
}

// Helper function to compute XC energy density
#[inline]
pub fn build_xc_density(
    name: impl Into<String>,
    density_matrix: Arc<dyn Expr>,
    overlap_distribution: Arc<dyn Expr>,
    order: u32,
) -> Result<Arc<dyn Expr>, TinnedError> {
    // Generalized density vector
    let density_vector = Trace::new(MatrixMul::new(vec![overlap_distribution, density_matrix])?)?;

    Ok(Composition::new(name, order, density_vector))
}

// Helper function to stream formatting to XC energy at a grid point
// Format: name[coef; factor1; factor2; ...], omit coef if one
#[inline]
pub fn fmt_xc_energy_term(f: &mut Formatter, name: &str, mul: &Mul) -> FmtResult {
    write!(f, "{}[", name)?;

    let mut first = true;

    let coef = mul.coefficient();
    if coef.is_one() {
        write!(f, "{}", coef)?;
        first = false;
    }

    for factor in mul.factors() {
        if !first {
            write!(f, "; ")?;
        }
        write!(f, "{}", factor)?;
        first = false;
    }

    write!(f, "]")
}

// Helper function to stream formatting to XC potential at a grid point
// Format: name[coef; factor1; factor2; ...], omit coef if one
#[inline]
pub fn fmt_xc_potential_term(f: &mut Formatter, name: &str, mul: &MatrixMul) -> FmtResult {
    write!(f, "{}[", name)?;

    let mut first = true;

    let coef = mul.coefficient();
    if !is_one_expr(coef) {
        write!(f, "{}", coef)?;
        first = false;
    }

    for factor in mul.factors() {
        if !first {
            write!(f, "; ")?;
        }
        write!(f, "{}", factor)?;
        first = false;
    }

    write!(f, "]")
}
