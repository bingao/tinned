use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::core::Expr;
use tinned::expressions::Number;
use tinned::public::generic_error;

use crate::c_support::{with_downcast_expr, with_downcast_val};
use crate::core::{CComplex64, CRational64, ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};
use crate::public::NumberToleranceHandle;

#[ffi_export]
pub extern "C" fn tinned_number_from_i64(n: i64) -> Option<ExprBox> {
    let expr_arc: Arc<dyn Expr> = Number::from_i64(n);
    Some(ExprBox::new(ExprHandle::new(expr_arc)))
}

#[ffi_export]
pub extern "C" fn tinned_number_from_f64(f: f64) -> Option<ExprBox> {
    let expr_arc: Arc<dyn Expr> = Number::from_f64(f);
    Some(ExprBox::new(ExprHandle::new(expr_arc)))
}

#[ffi_export]
pub extern "C" fn tinned_number_from_complex(z: CComplex64) -> Option<ExprBox> {
    let expr_arc: Arc<dyn Expr> = Number::from_complex(z.to_complex64());
    Some(ExprBox::new(ExprHandle::new(expr_arc)))
}

#[ffi_export]
pub extern "C" fn tinned_number_from_rational(
    q: CRational64,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let r = match q.to_rational64() {
        Some(v) => v,
        None => {
            tinned_error_new(
                out_err,
                generic_error("Denominator is zero in tinned_number_from_rational", None),
            );
            return None;
        },
    };

    let expr_arc: Arc<dyn Expr> = Number::from_rational(r);
    Some(ExprBox::new(ExprHandle::new(expr_arc)))
}

#[ffi_export]
pub extern "C" fn tinned_number_zero() -> Option<ExprBox> {
    let expr_arc: Arc<dyn Expr> = Number::zero();
    Some(ExprBox::new(ExprHandle::new(expr_arc)))
}

#[ffi_export]
pub extern "C" fn tinned_number_one() -> Option<ExprBox> {
    let expr_arc: Arc<dyn Expr> = Number::one();
    Some(ExprBox::new(ExprHandle::new(expr_arc)))
}

#[ffi_export]
pub extern "C" fn tinned_number_minus_one() -> Option<ExprBox> {
    let expr_arc: Arc<dyn Expr> = Number::minus_one();
    Some(ExprBox::new(ExprHandle::new(expr_arc)))
}

#[ffi_export]
pub extern "C" fn tinned_number_one_half() -> Option<ExprBox> {
    let expr_arc: Arc<dyn Expr> = Number::one_half();
    Some(ExprBox::new(ExprHandle::new(expr_arc)))
}

#[ffi_export]
pub extern "C" fn tinned_number_imaginary_unit() -> Option<ExprBox> {
    let expr_arc: Arc<dyn Expr> = Number::imaginary_unit();
    Some(ExprBox::new(ExprHandle::new(expr_arc)))
}

#[ffi_export]
pub extern "C" fn tinned_number_is_zero(
    h: Option<&ExprHandle>,
    tol: Option<&NumberToleranceHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    let num_tol = tol.map(NumberToleranceHandle::as_ref).cloned();
    with_downcast_val::<Number, bool>(h, out_err, "tinned_number_is_zero", |num| {
        num.is_zero(num_tol)
    })
    .unwrap_or(false)
}

#[ffi_export]
pub extern "C" fn tinned_number_is_one(
    h: Option<&ExprHandle>,
    tol: Option<&NumberToleranceHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    let num_tol = tol.map(NumberToleranceHandle::as_ref).cloned();
    with_downcast_val::<Number, bool>(h, out_err, "tinned_number_is_one", |num| num.is_one(num_tol))
        .unwrap_or(false)
}

impl_expr_getters!(
    Number;
    tinned_number_conjugate => |n| Ok(n.conjugate().into()),
    tinned_number_negate => |n| Ok(n.negate().into()),
);
