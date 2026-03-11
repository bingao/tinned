use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::TwoElecEnergy;
use tinned::perturbations::PertMultichain;
use tinned::public::generic_error;

use crate::c_support::tinned_string_from_cstr;
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::PertMultichainHandle;

#[ffi_export]
pub extern "C" fn tinned_two_elec_energy_new(
    name: Option<char_p::Ref<'_>>,
    inner_density: Option<&ExprHandle>,
    outer_density: Option<&ExprHandle>,
    allow_density_swap: bool,
    dependencies: Option<&PertMultichainHandle>,
    derivative: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(name) = tinned_string_from_cstr(name) else {
        tinned_error_new(
            out_err,
            generic_error("Null or invalid name passed to tinned_two_elec_energy_new", None),
        );
        return None;
    };

    let Some(inner_density) = inner_density else {
        tinned_error_new(
            out_err,
            generic_error("Null inner density passed to tinned_two_elec_energy_new", None),
        );
        return None;
    };
    let inner_density_arc = inner_density.clone_arc();

    let mut builder =
        TwoElecEnergy::builder(name, inner_density_arc).allow_density_swap(allow_density_swap);
    if let Some(outer_density) = outer_density {
        let outer_density_arc = outer_density.clone_arc();
        builder = builder.outer_density(outer_density_arc);
    }
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
    tinned_two_elec_energy_name : TwoElecEnergy => |op| op.name().to_string()
);

impl_expr_getters!(
    TwoElecEnergy;
    tinned_two_elec_energy_inner_density => |op| Ok(Arc::clone(op.inner_density())),
    tinned_two_elec_energy_outer_density => |op| Ok(Arc::clone(op.outer_density())),
);

impl_val_getters!(
    TwoElecEnergy;
    tinned_two_elec_energy_allow_density_swap: bool => |op| op.allow_density_swap(); default = false,
);

// Get `derivative` (cloned).
impl_pert_multichain_getter!(
    tinned_two_elec_energy_derivative : TwoElecEnergy => |op| op.derivative().clone()
);

// Get `dependencies` (cloned).
impl_pert_multichain_getter!(
    tinned_two_elec_energy_dependencies : TwoElecEnergy => |op| op.dependencies().clone()
);
