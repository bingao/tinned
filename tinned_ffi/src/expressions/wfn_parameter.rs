use safer_ffi::prelude::*;

use tinned::expressions::WfnParameter;
use tinned::public::generic_error;

use crate::c_support::tinned_string_from_cstr;
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};

#[ffi_export]
pub extern "C" fn tinned_wfn_parameter_new(
    name: Option<char_p::Ref<'_>>,
    is_perturbing: bool,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(name) = tinned_string_from_cstr(name) else {
        tinned_error_new(
            out_err,
            generic_error("Null or invalid name passed to tinned_wfn_parameter_new", None),
        );
        return None;
    };

    match <WfnParameter>::builder(name).is_perturbing(is_perturbing).build() {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Get `name` (caller must free the returned C string).
impl_cstr_getter!(
    tinned_wfn_parameter_name: WfnParameter => |obj| obj.name().to_string()
);

impl_val_getters!(
    WfnParameter;
    tinned_wfn_parameter_is_perturbing: bool => |obj| obj.is_perturbing(); default = false,
);

// Get `derivative` (cloned).
impl_pert_multichain_getter!(WfnParameter, tinned_wfn_parameter_derivative, |obj| obj
    .derivative()
    .clone());
