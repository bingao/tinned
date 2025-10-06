use std::{
    ffi::{CStr, CString},
    os::raw::c_char,
};

// Frees a `char*` that was allocated by Rust (via `CString::into_raw`). Safe
// to call with `NULL`.
#[unsafe(no_mangle)]
pub extern "C" fn tinned_free_cstring(s: *mut c_char) {
    if s.is_null() {
        return;
    }
    unsafe {
        let _ = CString::from_raw(s);
    }
}

// Converts a `NULL`-terminated const `char*` from C into a Rust String.
// Returns `None` if `ptr` is `NULL`.
pub(crate) fn cstr_to_string(ptr: *const c_char) -> Option<String> {
    if ptr.is_null() {
        return None;
    }
    // SAFETY: `ptr` must point to a valid, `NULL`-terminated C string.
    // The caller (FFI boundary) is responsible for ensuring this.
    unsafe { Some(CStr::from_ptr(ptr).to_string_lossy().into_owned()) }
}

// Convert a Rust string (`&str` or `String`) into a heap-allocated,
// `NULL`-terminated `char*` that C can read. C must free with
// `tinned_free_cstring()`.
pub(crate) fn to_cstring<S: AsRef<str>>(s: S) -> *mut c_char {
    CString::new(s.as_ref()).map(CString::into_raw).unwrap_or(std::ptr::null_mut())
}
