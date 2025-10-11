use safer_ffi::prelude::*;

use tinned::core::TinnedError;

/// An *opaque* handle that C can only pass around
#[derive_ReprC]
#[repr(opaque)]
pub struct TinnedErrorHandle {
    inner: TinnedError,
}

/// Owned box by C after return
pub type TinnedErrorBox = repr_c::Box<TinnedErrorHandle>;

impl TinnedErrorHandle {
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
}

// Internal: fill an out-err safer-ffi out-parameter (used by FFI entrypoints)
#[inline]
pub(crate) fn set_out_err(out_err: Option<Out<'_, TinnedErrorBox>>, err: TinnedError) {
    if let Some(out) = out_err {
        let rust_box = std::boxed::Box::new(TinnedErrorHandle::new(err));
        out.write(rust_box.into());
    }
}
