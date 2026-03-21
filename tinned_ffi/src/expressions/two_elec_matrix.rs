use safer_ffi::prelude::*;
use std::collections::BTreeSet;
use std::sync::Arc;

use tinned::expressions::TwoElecMatrix;
use tinned::perturbations::Perturbation;
use tinned::public::generic_error;

use crate::c_support::{ffi_map_expr_as_pertvec, tinned_string_from_cstr};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::{
    PertMultichainHandle, PerturbationBox, PerturbationSlice, perturbation_set_from_slice,
};

#[ffi_export]
pub extern "C" fn tinned_two_elec_matrix_new(
    name: Option<char_p::Ref<'_>>,
    is_perturbing: bool,
    dependencies: Option<&PertMultichainHandle>,
    independent_perturbations: Option<&PerturbationSlice>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(name) = tinned_string_from_cstr(name) else {
        tinned_error_new(
            out_err,
            generic_error("Null or invalid name passed to tinned_two_elec_matrix_new", None),
        );
        return None;
    };

    let mut builder = TwoElecMatrix::builder(name).is_perturbing(is_perturbing);

    if let Some(dependencies) = dependencies {
        let deps = dependencies.as_ref().clone();
        builder = builder.dependencies(deps);
    }

    if let Some(slice) = independent_perturbations {
        let indep_perts = match perturbation_set_from_slice::<BTreeSet<Arc<Perturbation>>>(
            slice,
            "tinned_two_elec_matrix_new",
        ) {
            Ok(set) => set,
            Err(err) => {
                tinned_error_new(out_err, err);
                return None;
            },
        };

        builder = builder.independent_perturbations(indep_perts);
    }

    match builder.build() {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(err) => {
            tinned_error_new(out_err, err);
            None
        },
    }
}

// Get `name` (caller must free the returned C string).
impl_cstr_getter!(
    tinned_two_elec_matrix_name: TwoElecMatrix => |op| op.name().to_string()
);

impl_val_getters!(
    TwoElecMatrix;
    tinned_two_elec_matrix_is_perturbing: bool => |op| op.is_perturbing(); default = false,
);

// Get `dependencies` (cloned).
impl_pert_multichain_getter!(TwoElecMatrix, tinned_two_elec_matrix_dependencies, |op| op
    .dependencies()
    .clone());

// Get independent perturbations as `repr_c::Vec<PerturbationBox>`, must be
// freed by calling `tinned_perturbation_vec_free()`
impl_vec_getter!(
    TwoElecMatrix,
    tinned_two_elec_matrix_independent_perturbations,
    repr_c::Vec<PerturbationBox>,
    ffi_map_expr_as_pertvec,
    |op| op.independent_perturbations()
);

// Get `derivative` (cloned).
impl_pert_multichain_getter!(TwoElecMatrix, tinned_two_elec_matrix_derivative, |op| op
    .derivative()
    .clone());
