use std::error::Error as StdError;

use tinned::core::TinnedError;

use crate::c_support::to_cstring;
use crate::core::TinnedErrorBox;

// Public: free the error box
#[unsafe(no_mangle)]
pub extern "C" fn tinned_error_free(err: *mut TinnedErrorBox) {
    if err.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(err));
    }
}

// Build a human-readable string:
// - First line: error message
// - Then the source() chain, each on its own line.
fn format_error(e: &TinnedError) -> String {
    use std::fmt::Write;

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

// Get printable error text. C must free with `tinned_free_cstring()`.
#[unsafe(no_mangle)]
pub extern "C" fn tinned_error_display(err: *const TinnedErrorBox) -> *mut std::os::raw::c_char {
    if err.is_null() {
        return to_cstring("No error");
    }
    let e = unsafe { &*err }.as_ref();
    to_cstring(format_error(e))
}
