use safer_ffi::prelude::*;

use tinned::expressions::OneElecOperator;
use tinned::perturbations::PertMultichain;
use tinned::public::generic_error;

use crate::c_support::tinned_string_from_cstr;
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::perturbations::PertMultichainHandle;

#[ffi_export]
pub extern "C" fn tinned_one_elec_operator_new(
    name: Option<char_p::Ref<'_>>,
    dependencies: Option<&PertMultichainHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let Some(name) = tinned_string_from_cstr(name) else {
        tinned_error_new(
            out_err,
            generic_error("Null or invalid name passed to tinned_one_elec_operator_new", None),
        );
        return None;
    };

    let mut builder = OneElecOperator::builder(name);
    if let Some(dependencies) = dependencies {
        let deps: PertMultichain = dependencies.as_ref().clone();
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
    tinned_one_elec_operator_name : OneElecOperator => |op| op.name().to_string()
);

// Get `derivative` (cloned).
impl_pert_multichain_getter!(
    tinned_one_elec_operator_derivative : OneElecOperator => |op| op.derivative().clone()
);

// Get `dependencies` (cloned).
impl_pert_multichain_getter!(
    tinned_one_elec_operator_dependencies : OneElecOperator => |op| op.dependencies().clone()
);
