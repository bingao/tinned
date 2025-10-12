use safer_ffi::prelude::*;
use std::error::Error as StdError;
use std::fmt::Write;

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
    pub fn new(err: TinnedError) -> Self {
        Self {
            inner: err,
        }
    }

    #[inline]
    pub fn as_ref(&self) -> &TinnedError {
        &self.inner
    }
}

// Creates an opaque handle for an error and fill into a box (used by FFI entrypoints).
#[inline]
pub fn tinned_error_new(out_err: Option<Out<'_, TinnedErrorBox>>, err: TinnedError) {
    if let Some(out) = out_err {
        out.write(TinnedErrorBox::new(TinnedErrorHandle::new(err)));
    }
}

// Public: free the error box
#[ffi_export]
pub fn tinned_error_free(err: Option<TinnedErrorBox>) {
    drop(err);
}

// Build a human-readable string:
// - First line: error message
// - Then the source() chain, each on its own line.
fn format_error(e: &TinnedError) -> String {
    let mut s = e.to_string();
    let mut cur: Option<&dyn StdError> = (e as &dyn StdError).source();
    let mut idx = 1;
    if cur.is_some() {
        s.push_str("\nCaused by:");
    }
    while let Some(src) = cur {
        let _ = write!(&mut s, "\n  [{}] {}", idx, src);
        cur = src.source();
        idx += 1;
    }
    s
}

// Get printable error text. C must free with `tinned_string_free()`.
#[ffi_export]
pub fn tinned_error_display(err: Option<&TinnedErrorHandle>) -> char_p::Box {
    match err {
        Some(h) => char_p::new(format_error(h.as_ref())),
        None => char_p::new(String::new()),
    }
}
