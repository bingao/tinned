use std::{os::raw::c_char, ptr::null_mut, sync::Arc};

use tinned::core::{Expr, TinnedError};
use tinned::public::expression_error;

use crate::c_support::to_cstring;
use crate::core::{ExprBox, TinnedErrorBox, set_out_err, with_expr_or_err};

#[inline]
fn invalid_type_err<T: 'static>(caller: &'static str, expr: &Arc<dyn Expr>) -> TinnedError {
    let msg: &'static str = Box::leak(
        format!("Invalid expression passed to {caller}; expected {}", std::any::type_name::<T>())
            .into_boxed_str(),
    );
    expression_error(msg, expr, None)
}

// (1) Primitive / copy types; use Default as the fallback to avoid threading a literal.
#[inline]
pub(crate) fn with_downcast_val<T: 'static, R: Copy + Default>(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
    caller: &'static str,
    f: impl FnOnce(&T) -> R,
) -> R {
    with_expr_or_err(h, out_err, caller, |e: &Arc<dyn Expr>| {
        if let Some(t) = e.as_any().downcast_ref::<T>() {
            f(t)
        } else {
            set_out_err(out_err, invalid_type_err::<T>(caller, e));
            R::default()
        }
    })
    .unwrap_or_default()
}

// (2) C string
#[inline]
pub(crate) fn with_downcast_cstr<T: 'static>(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
    caller: &'static str,
    f: impl FnOnce(&T) -> String,
) -> *mut c_char {
    with_expr_or_err(h, out_err, caller, |e: &Arc<dyn Expr>| {
        if let Some(t) = e.as_any().downcast_ref::<T>() {
            to_cstring(f(t))
        } else {
            set_out_err(out_err, invalid_type_err::<T>(caller, e));
            null_mut()
        }
    })
    .unwrap_or(null_mut())
}

// (3) Always accept a Result; For no-error paths, just return Ok(...)
#[inline]
pub(crate) fn with_downcast_expr_res<T: 'static>(
    h: *const ExprBox,
    out_err: *mut *mut TinnedErrorBox,
    caller: &'static str,
    f: impl FnOnce(&T) -> Result<Arc<dyn Expr>, TinnedError>,
) -> *mut ExprBox {
    with_expr_or_err(h, out_err, caller, |e: &Arc<dyn Expr>| {
        if let Some(t) = e.as_any().downcast_ref::<T>() {
            match f(t) {
                Ok(expr) => ExprBox::new(expr).into_raw(),
                Err(err) => {
                    set_out_err(out_err, err);
                    null_mut()
                },
            }
        } else {
            set_out_err(out_err, invalid_type_err::<T>(caller, e));
            null_mut()
        }
    })
    .unwrap_or(null_mut())
}
