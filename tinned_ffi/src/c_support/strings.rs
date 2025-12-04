use safer_ffi::prelude::*;

// Unified free for strings returned from FFI (`char_p::Box`)
#[cfg(feature = "tinned-ffi")]
#[ffi_export]
pub fn tinned_string_free(s: Option<char_p::Box>) {
    drop(s);
}

// C `const char*` -> `Option<String>``. `None`` if `NULL` or invalid UTF-8.
#[inline]
pub(crate) fn tinned_string_from_cstr(s: Option<char_p::Ref<'_>>) -> Option<String> {
    s.map(|r| r.to_str().to_owned())
}

// Rust string -> owned `char*` for C (free with `tinned_string_free`).
#[inline]
pub(crate) fn tinned_string_to_cstr<S: AsRef<str>>(s: S) -> char_p::Box {
    char_p::new(s.as_ref())
}
