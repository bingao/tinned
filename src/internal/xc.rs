use std::sync::Arc;

use crate::core::{Expr, TinnedError};
use crate::expressions::{
    Composition, MatrixMul, NonElecFunction, OneElecOperator, Trace, WfnParameter,
};
use crate::public::{expression_error, is_expr_type};

// Helper function to validate the density matrix, grid weight and overlap
// distribution
#[inline]
pub(crate) fn validate_xc_inputs(
    density_matrix: &Arc<dyn Expr>,
    grid_weight: &Arc<dyn Expr>,
    overlap_distribution: &Arc<dyn Expr>,
) -> Result<(), TinnedError> {
    if !is_expr_type::<WfnParameter>(density_matrix) {
        return Err(expression_error(
            "validate_xc_inputs() - density matrix must be WfnParameter",
            density_matrix,
            None,
        ));
    }

    if !is_expr_type::<NonElecFunction>(grid_weight) {
        return Err(expression_error(
            "validate_xc_inputs() - grid weight must be NonElecFunction",
            grid_weight,
            None,
        ));
    }

    if !is_expr_type::<OneElecOperator>(overlap_distribution) {
        return Err(expression_error(
            "validate_xc_inputs() - overlap distribution must be OneElecOperator",
            overlap_distribution,
            None,
        ));
    }

    Ok(())
}

// Helper function to compute XC energy density
#[inline]
pub(crate) fn build_xc_density(
    name: impl Into<String>,
    density_matrix: Arc<dyn Expr>,
    overlap_distribution: Arc<dyn Expr>,
    order: u32,
) -> Result<Arc<dyn Expr>, TinnedError> {
    // Generalized density vector
    let density_vector = Trace::new(MatrixMul::new(vec![overlap_distribution, density_matrix])?)?;

    Composition::new(name, order, density_vector)
}
