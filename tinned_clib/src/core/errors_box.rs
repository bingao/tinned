use tinned::core::TinnedError;

// Opaque error box for FFI
pub struct TinnedErrorBox {
    inner: TinnedError,
}

impl TinnedErrorBox {
    #[inline]
    pub(crate) fn new(err: TinnedError) -> Self {
        Self {
            inner: err,
        }
    }

    #[inline]
    pub(crate) fn as_ref(&self) -> &TinnedError {
        &self.inner
    }

    // Turn this box into a raw pointer for FFI returns
    #[inline]
    pub(crate) fn into_raw(self) -> *mut TinnedErrorBox {
        Box::into_raw(Box::new(self))
    }
}

// Internal: fill an out-err pointer (used by FFI entrypoints)
#[inline]
pub(crate) fn set_out_err(out_err: *mut *mut TinnedErrorBox, err: TinnedError) {
    if out_err.is_null() {
        return;
    }
    let h = TinnedErrorBox::new(err).into_raw();
    // With `#![deny(unsafe_op_in_unsafe_fn)]`, keep an explicit unsafe block.
    unsafe {
        *out_err = h;
    }
}
