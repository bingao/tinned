use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::expressions::WfnParameter;
use tinned::public::generic_error;

use crate::c_support::{ffi_map_expr_as, tinned_string_from_cstr, tinned_string_to_cstr};
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::{PertMultichainBox, PertMultichainHandle};

#[ffi_export]
pub extern "C" fn tinned_wfn_parameter_new(
    name: Option<char_p::Ref<'_>>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(name) = tinned_string_from_cstr(name) else {
        tinned_error_new(
            out_err,
            generic_error("Null or invalid name passed to tinned_wfn_parameter_new", None),
        );
        return None;
    };

    match <WfnParameter>::builder(name).build() {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Get `name` (caller must free the returned C string).
impl_cstr_getter!(
    tinned_wfn_parameter_name : WfnParameter => |lag| lag.name().to_string()
);

// Get `derivative` (cloned).
impl_pert_multichain_getter!(
    tinned_wfn_parameter_derivative : WfnParameter => |obj| obj.derivative().clone()
);
