use safer_ffi::prelude::*;
use std::error::Error as StdError;
use std::fmt::Write;

use tinned::core::TinnedError;

use crate::core::{TinnedErrorBox, TinnedErrorHandle};

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
