use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::NonElecFunction;
use tinned::perturbations::PertMultichain;
use tinned::public::generic_error;

use crate::c_support::{
    tinned_string_from_cstr, with_downcast_cstr, with_downcast_pert_multichain,
};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::{PertMultichainBox, PertMultichainHandle};

#[ffi_export]
pub extern "C" fn tinned_non_elec_function_new(
    name: Option<char_p::Ref<'_>>,
    dependencies: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(name) = tinned_string_from_cstr(name) else {
        tinned_error_new(
            out_err,
            generic_error("Null or invalid name passed to tinned_non_elec_function_new", None),
        );
        return None;
    };

    let mut builder = NonElecFunction::builder(name);
    if let Some(handle) = dependencies {
        let deps: PertMultichain = handle.as_ref().clone();
        builder = builder.dependencies(deps);
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
    tinned_non_elec_function_name : NonElecFunction => |op| op.name().to_string()
);

// Get `derivative` (cloned).
impl_pert_multichain_getter!(
    tinned_non_elec_function_derivative : NonElecFunction => |op| Ok(Arc::new(op.derivative().clone()))
);

// Get `dependencies` (cloned).
impl_pert_multichain_getter!(
    tinned_non_elec_function_dependencies : NonElecFunction => |op| Ok(Arc::new(op.dependencies().clone()))
);
