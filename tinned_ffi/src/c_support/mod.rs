mod expr_downcast;
mod pointers;
mod strings;

pub use pointers::{
    try_from_handle, try_map_from_slices, try_set_from_slice, try_vec_from_slice, try_with_handle,
};

pub(crate) use expr_downcast::{
    with_downcast_cstr, with_downcast_expr, with_downcast_pert, with_downcast_pert_multichain,
    with_downcast_val,
};
pub(crate) use strings::{tinned_string_from_cstr, tinned_string_to_cstr};
