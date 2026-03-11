use safer_ffi::prelude::*;

use tinned::expressions::Mul;

use crate::c_support::ffi_map_expr_as_exprvec;
use crate::core::{
    ExprBox, ExprHandle, ExprSlice, TinnedErrorBox, expr_vec_from_slice, tinned_error_new,
};

#[ffi_export]
pub extern "C" fn tinned_mul_new(
    terms: &ExprSlice,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let terms_vec = match expr_vec_from_slice(terms, "tinned_mul_new") {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            return None;
        },
    };

    match Mul::new(terms_vec) {
        Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

impl_expr_getters!(
    Mul;
    tinned_mul_coefficient => |mul| Ok(mul.coefficient().into()),
);

// Returns a cloned vector of factors
#[ffi_export]
pub extern "C" fn tinned_mul_factors(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> repr_c::Vec<ExprBox> {
    ffi_map_expr_as_exprvec::<Mul>(h, out_err, "tinned_mul_factors", |mul| mul.factors())
}
