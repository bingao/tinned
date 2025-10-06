use std::{ptr::null_mut, sync::Arc};

use tinned::expressions::Conjugate;
use tinned::public::generic_error;

use crate::c_support::with_downcast_expr_res;
use crate::core::{ExprBox, TinnedErrorBox, set_out_err};

#[unsafe(no_mangle)]
pub extern "C" fn tinned_conjugate_new(
    argument_ptr: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
    if argument_ptr.is_null() {
        set_out_err(out_err, generic_error("Null argument passed to tinned_conjugate_new", None));
        return null_mut();
    }
    let argument = unsafe { (&*argument_ptr).arc_clone() };

    match Conjugate::new(argument) {
        Ok(expr) => ExprBox::new(expr).into_raw(),
        Err(e) => {
            set_out_err(out_err, e);
            null_mut()
        },
    }
}

// Get `argument` (cloned).
#[unsafe(no_mangle)]
pub extern "C" fn tinned_conjugate_argument(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
) -> *mut ExprBox {
    with_downcast_expr_res::<Conjugate>(h, out_err, "tinned_conjugate_argument", |cj| {
        Ok(Arc::clone(cj.argument()))
    })
}
