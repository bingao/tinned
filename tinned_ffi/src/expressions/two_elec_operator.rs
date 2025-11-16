use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::TwoElecOperator;
use tinned::perturbations::PertMultichain;
use tinned::public::generic_error;

use crate::c_support::{
    tinned_string_from_cstr, ffi_map_expr_as, tinned_string_to_cstr,
};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::{PertMultichainBox, PertMultichainHandle};

#[ffi_export]
pub extern "C" fn tinned_two_elec_operator_new(
    name: Option<char_p::Ref<'_>>,
    density: Option<&ExprHandle>,
    dependencies: Option<&PertMultichainHandle>,
    derivative: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(name) = tinned_string_from_cstr(name) else {
        tinned_error_new(
            out_err,
            generic_error("Null or invalid name passed to tinned_two_elec_operator_new", None),
        );
        return None;
    };

    let Some(density) = density else {
        tinned_error_new(
            out_err,
            generic_error("Null density passed to tinned_two_elec_operator_new", None),
        );
        return None;
    };
    let density_arc = density.clone_arc();

    let mut builder = TwoElecOperator::builder(name, density_arc);
    if let Some(dependencies) = dependencies {
        let deps: PertMultichain = dependencies.as_ref().clone();
        builder = builder.dependencies(deps);
    }
    if let Some(derivative) = derivative {
        let deriv: PertMultichain = derivative.as_ref().clone();
        builder = builder.derivative(deriv);
    }

    match builder.build() {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Get `name` (caller must free the returned C string).
impl_cstr_getter!(
    tinned_two_elec_operator_name : TwoElecOperator => |op| op.name().to_string()
);

impl_expr_getters!(
    TwoElecOperator;
    tinned_two_elec_operator_density => |op| Ok(Arc::clone(op.density())),
);

// Get `derivative` (cloned).
impl_pert_multichain_getter!(
    tinned_two_elec_operator_derivative : TwoElecOperator => |op| op.derivative().clone()
);

// Get `dependencies` (cloned).
impl_pert_multichain_getter!(
    tinned_two_elec_operator_dependencies : TwoElecOperator => |op| op.dependencies().clone()
);
