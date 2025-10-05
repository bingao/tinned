use tinned::core::TinnedError;

// Opaque error box for FFI
pub struct TinnedErrorBox {
    inner: TinnedError,
}

impl TinnedErrorBox {
    #[inline]
    pub(crate) fn as_ref(&self) -> &TinnedError {
        &self.inner
    }
}

// Internal: box -> raw pointer
#[inline]
fn error_box_from(err: TinnedError) -> *mut TinnedErrorBox {
    Box::into_raw(Box::new(TinnedErrorBox {
        inner: err,
    }))
}

// Internal: fill an out-err pointer (used by FFI entrypoints)
#[inline]
pub(crate) fn set_out_err(out_err: *mut *mut TinnedErrorBox, err: TinnedError) {
    if out_err.is_null() {
        return;
    }
    unsafe {
        *out_err = error_box_from(err);
    }
}
