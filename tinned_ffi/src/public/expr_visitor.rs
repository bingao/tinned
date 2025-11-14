use safer_ffi::prelude::*;
use std::sync::Arc;

use tinned::core::{Expr, TinnedError};
use tinned::public::{ExprTag, ExprVisitor, generic_error, walk_expr_postorder};

use crate::c_support::try_from_handle;
use crate::core::{ExprBox, ExprHandle, TinnedErrorBox, tinned_error_new};

#[repr(C)]
#[derive_ReprC]
pub struct CExprVisitor {
    pub ctx: *mut core::ffi::c_void,
    pub begin_node: extern "C" fn(*mut core::ffi::c_void, ExprTag, usize) -> bool,
    // Owned handle transferred to C; C must free it with tinned_expr_free.
    pub on_leaf: extern "C" fn(*mut core::ffi::c_void, ExprTag, ExprBox) -> bool,
    pub end_node: extern "C" fn(*mut core::ffi::c_void, ExprTag, usize) -> bool,
}

struct CVisitorBridge {
    visitor: CExprVisitor,
}

impl ExprVisitor for CVisitorBridge {
    fn begin(&mut self, tag: ExprTag, arity: usize) -> Result<(), TinnedError> {
        if (self.visitor.begin_node)(self.visitor.ctx, tag, arity) {
            Ok(())
        } else {
            Err(generic_error("begin_node returned false", None))
        }
    }

    fn leaf(&mut self, tag: ExprTag, expr: &Arc<dyn Expr>) -> Result<(), TinnedError> {
        // Build an owned handle for C. C must free it.
        let h = ExprBox::new(ExprHandle::new(Arc::clone(expr)));
        if (self.visitor.on_leaf)(self.visitor.ctx, tag, h) {
            Ok(())
        } else {
            Err(generic_error("on_leaf returned false", None))
        }
    }

    fn end(&mut self, tag: ExprTag, arity: usize) -> Result<(), TinnedError> {
        if (self.visitor.end_node)(self.visitor.ctx, tag, arity) {
            Ok(())
        } else {
            Err(generic_error("end_node returned false", None))
        }
    }
}

#[ffi_export]
pub extern "C" fn tinned_walk_expr_postorder(
    h: Option<&ExprHandle>,
    visitor: CExprVisitor,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    let expr =
        match try_from_handle(h, "tinned_expr_eval_with_c", "ExprHandle", |eh| eh.clone_arc()) {
            Ok(x) => x,
            Err(e) => {
                tinned_error_new(out_err, e);
                return false;
            },
        };
    let mut bridge = CVisitorBridge {
        visitor,
    };
    match walk_expr_postorder(&expr, &mut bridge) {
        Ok(()) => true,
        Err(e) => {
            tinned_error_new(out_err, e);
            false
        },
    }
}
