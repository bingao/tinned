mod expr_downcast;
mod pointers;
mod strings;

pub use pointers::with_box_or_err;

pub(crate) use expr_downcast::{with_downcast_cstr, with_downcast_expr_res, with_downcast_val};
pub(crate) use pointers::vec_arc_from_ptrs;
pub(crate) use strings::{cstr_to_string, to_cstring};
