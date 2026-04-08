use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::BasisTimeEvolution;
use tinned::perturbations::PertMultichain;
use tinned::public::generic_error;

use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::PertMultichainHandle;

#[ffi_export]
pub extern "C" fn tinned_basis_time_evolution_new(
    dependencies: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(dependencies) = dependencies else {
        tinned_error_new(
            out_err,
            generic_error("Null dependencies passed to tinned_basis_time_evolution_new", None),
        );
        return None;
    };
    let deps: PertMultichain = dependencies.as_ref().clone();

    match BasisTimeEvolution::builder(deps).build() {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

impl_val_getters!(
    BasisTimeEvolution;
    tinned_basis_time_evolution_at_zero_perturbations: bool => |op| op.at_zero_perturbations(); default = false,
);

impl_expr_getters!(
    BasisTimeEvolution;
    tinned_basis_time_evolution_braket => |op| Ok(Arc::clone(op.braket())),
);

// Get `derivative` (cloned).
impl_pert_multichain_getter!(BasisTimeEvolution, tinned_basis_time_evolution_derivative, |op| op
    .derivative()
    .clone());

// Get `dependencies` (cloned).
impl_pert_multichain_getter!(BasisTimeEvolution, tinned_basis_time_evolution_dependencies, |op| op
    .dependencies()
    .clone());
