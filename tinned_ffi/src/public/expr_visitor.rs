#[repr(C)]
#[derive(safer_ffi::derive_ReprC)]
pub struct CExprVisitor {
    pub ctx: *mut core::ffi::c_void,
    pub begin_node: extern "C" fn(*mut core::ffi::c_void, ExprTag, usize) -> bool,
    // Owned handle transferred to C; C must free it.
    pub on_leaf: extern "C" fn(*mut core::ffi::c_void, ExprTag, ExprHandle) -> bool,
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
            Err(TinnedError::msg_static("begin_node returned false"))
        }
    }

    fn leaf(&mut self, tag: ExprTag, expr: &Arc<dyn Expr>) -> Result<(), TinnedError> {
        // Build an owned handle for C. C must free it.
        let h = ExprHandle::new(Arc::clone(expr));
        if (self.visitor.on_leaf)(self.visitor.ctx, tag, h) {
            Ok(())
        } else {
            Err(TinnedError::new_static("on_leaf returned false"))
        }
    }

    fn end(&mut self, tag: ExprTag, arity: usize) -> Result<(), TinnedError> {
        if (self.visitor.end_node)(self.visitor.ctx, tag, arity) {
            Ok(())
        } else {
            Err(TinnedError::msg_static("end_node returned false"))
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
