macro_rules! impl_cstr_getter {
    ($fn_name:ident : $type_name:path => |$obj:ident| $body:expr) => {
        #[ffi_export]
        pub extern "C" fn $fn_name(
            h: Option<&ExprHandle>,
            out_err: Option<Out<'_, TinnedErrorBox>>,
        ) -> Option<char_p::Box> {
            ffi_map_expr_as::<$type_name, _>(h, out_err, stringify!($fn_name), |$obj| {
                Ok(tinned_string_to_cstr($body))
            })
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
                ffi_map_expr_as::<$type_name, _>(
                    h,
                    out_err,
                    stringify!($fn_name),
                    |$obj| {
                        $body.map(|arc| ExprBox::new(ExprHandle::new(arc)))
                    }
                )
            }
        )+
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
                ffi_map_expr_as_copy::<$type_name, $return_type>(
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
            ffi_map_expr_as::<$type_name, _>(
                h,
                out_err,
                stringify!($fn_name),
                |$obj| {
                    let chain = Arc::new($body);
                    Ok(PertMultichainBox::new(PertMultichainHandle::new(chain)))
                },
            )

        }
    };
}
