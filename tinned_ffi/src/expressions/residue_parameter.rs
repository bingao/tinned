use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::ResidueParameter;
use tinned::public::generic_error;

use crate::PerturbationHandle;
use crate::c_support::{ffi_map_expr_as, ffi_map_expr_as_copy};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::{PerturbationBox, PerturbationSlice, perturbation_vec_from_slice};

#[ffi_export]
pub extern "C" fn tinned_residue_parameter_new(
    perturbations: PerturbationSlice<'_>,
    excited_state: Option<&ExprHandle>,
    parameter: Option<&ExprHandle>,
    positive_frequency: bool,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let vec = match perturbation_vec_from_slice(perturbations, "tinned_residue_parameter_new") {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            return None;
        },
    };

    let Some(excited_state) = excited_state else {
        tinned_error_new(
            out_err,
            generic_error("Null excited state passed to tinned_residue_parameter_new", None),
        );
        return None;
    };
    let excited_state_arc = excited_state.clone_arc();

    let Some(parameter) = parameter else {
        tinned_error_new(
            out_err,
            generic_error("Null parameter passed to tinned_residue_parameter_new", None),
        );
        return None;
    };
    let parameter_arc = parameter.clone_arc();

    match ResidueParameter::builder(vec, excited_state_arc, parameter_arc)
        .positive_frequency(positive_frequency)
        .build()
    {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

impl_val_getters!(
    ResidueParameter;
    tinned_residue_parameter_positive_frequency: bool => |res| res.positive_frequency(); default = false,
);

// Returns a cloned vector of perturbations
#[ffi_export]
pub extern "C" fn tinned_residue_parameter_perturbations(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> repr_c::Vec<PerturbationBox> {
    match ffi_map_expr_as::<ResidueParameter, _>(h, out_err, "tinned_residue_parameter_perturbations", |res| {
        let mut perturbations: Vec<PerturbationBox> = Vec::with_capacity(res.perturbations().len());
        for pert in res.perturbations() {
            perturbations.push(PerturbationBox::new(PerturbationHandle::new(Arc::clone(pert))));
        }
        Ok(perturbations.into())
    }) {
        Some(perturbations) => perturbations,
        None => Vec::<PerturbationBox>::new().into(),
    }
}

impl_expr_getters!(
    ResidueParameter;
    tinned_residue_parameter_excited_state => |res| Ok(Arc::clone(res.excited_state())),
    tinned_residue_parameter_parameter => |res| Ok(Arc::clone(res.parameter())),
);
