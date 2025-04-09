use std::sync::Arc;

use crate::core::{Expr, TinnedError};
use crate::expressions::{Composition, MatrixMul, Trace, WfnParameter};
use crate::utils::invalid_expression_error;

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
