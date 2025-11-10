use safer_ffi::prelude::*;

use tinned::expressions::ZeroOperator;

use crate::core::{ExprBox, ExprHandle};

/// Create a new `ZeroOperator` expression.
#[ffi_export]
pub extern "C" fn tinned_zero_operator_new() -> Option<ExprBox> {
    Some(ExprBox::new(ExprHandle::new(ZeroOperator::new())))
}
