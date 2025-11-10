macro_rules! impl_cstr_getter {
    ($fn_name:ident : $type_name:path => |$obj:ident| $body:expr) => {
        #[ffi_export]
        pub extern "C" fn $fn_name(
            h: Option<&ExprHandle>,
            out_err: Option<Out<'_, TinnedErrorBox>>,
        ) -> Option<char_p::Box> {
            with_downcast_cstr::<$type_name>(h, out_err, stringify!($fn_name), |$obj| $body)
        }
    };
}

macro_rules! impl_expr_getters {
    ($type_name:path; $( $fn_name:ident => |$obj:ident| $body:expr ),+ $(,)?) => {
        $(
            #[ffi_export]
            pub extern "C" fn $fn_name(
                h: Option<&ExprHandle>,
                out_err: Option<Out<'_, TinnedErrorBox>>,
            ) -> Option<ExprBox> {
                with_downcast_expr::<$type_name>(
                    h,
                    out_err,
                    stringify!($fn_name),
                    |$obj| $body
                )
            }
        )+
    };
}

macro_rules! impl_expr_index_getter {
    ($fn_name:ident : $type_name:path => $collection:ident) => {
        #[ffi_export]
        pub extern "C" fn $fn_name(
            h: Option<&ExprHandle>,
            i: usize,
            out_err: Option<Out<'_, TinnedErrorBox>>,
        ) -> Option<ExprBox> {
            with_downcast_expr::<$type_name>(h, out_err, stringify!($fn_name), |obj| {
                obj.$collection().get(i).cloned().ok_or_else(|| {
                    generic_error(
                        format!(
                            "Index {} out of bounds (len = {}) in {}",
                            i,
                            obj.$collection().len(),
                            stringify!($fn_name)
                        ),
                        None,
                    )
                })
            })
        }
    };
}

macro_rules! impl_val_getters {
    ($type_name:path; $( $fn_name:ident : $return_type:ty => |$obj:ident| $body:expr ; default = $def:expr ),+ $(,)?) => {
        $(
            #[ffi_export]
            pub extern "C" fn $fn_name(
                h: Option<&ExprHandle>,
                out_err: Option<Out<'_, TinnedErrorBox>>,
            ) -> $return_type {
                with_downcast_val::<$type_name, $return_type>(
                    h,
                    out_err,
                    stringify!($fn_name),
                    |$obj| { $body }
                )
                .unwrap_or($def)
            }
        )+
    };
}

macro_rules! impl_pert_multichain_getter {
    ($fn_name:ident : $type_name:path => |$obj:ident| $body:expr) => {
        #[ffi_export]
        pub extern "C" fn $fn_name(
            h: Option<&ExprHandle>,
            out_err: Option<Out<'_, TinnedErrorBox>>,
        ) -> Option<PertMultichainBox> {
            with_downcast_pert_multichain::<$type_name>(h, out_err, stringify!($fn_name), |$obj| {
                $body
            })
        }
    };
}
