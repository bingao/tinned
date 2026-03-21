macro_rules! impl_cstr_getter {
    ($fn_name:ident : $type_name:path => |$obj:ident| $body:expr) => {
        #[::safer_ffi::ffi_export]
        pub extern "C" fn $fn_name(
            h: ::std::option::Option<&$crate::core::ExprHandle>,
            out_err: ::std::option::Option<
                ::safer_ffi::prelude::Out<'_, $crate::core::TinnedErrorBox>,
            >,
        ) -> ::std::option::Option<::safer_ffi::prelude::char_p::Box> {
            $crate::c_support::ffi_map_expr_as::<$type_name, _>(
                h,
                out_err,
                stringify!($fn_name),
                |$obj| Ok($crate::c_support::tinned_string_to_cstr($body)),
            )
        }
    };
}

macro_rules! impl_expr_getters {
    ($type_name:path; $( $fn_name:ident => |$obj:ident| $body:expr ),+ $(,)?) => {
        $(
            #[::safer_ffi::ffi_export]
            pub extern "C" fn $fn_name(
                h: ::std::option::Option<&$crate::core::ExprHandle>,
                out_err: ::std::option::Option<::safer_ffi::prelude::Out<'_, $crate::core::TinnedErrorBox>>,
            ) -> ::std::option::Option<$crate::core::ExprBox> {
                $crate::c_support::ffi_map_expr_as::<$type_name, _>(
                    h,
                    out_err,
                    stringify!($fn_name),
                    |$obj| {
                        $body.map(|arc| $crate::core::ExprBox::new($crate::core::ExprHandle::new(arc)))
                    }
                )
            }
        )+
    };
}

macro_rules! impl_val_getters {
    ($type_name:path; $( $fn_name:ident : $return_type:ty => |$obj:ident| $body:expr ; default = $def:expr ),+ $(,)?) => {
        $(
            #[::safer_ffi::ffi_export]
            pub extern "C" fn $fn_name(
                h: ::std::option::Option<&$crate::core::ExprHandle>,
                out_err: ::std::option::Option<::safer_ffi::prelude::Out<'_, $crate::core::TinnedErrorBox>>,
            ) -> $return_type {
                $crate::c_support::ffi_map_expr_as_copy::<$type_name, $return_type>(
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

macro_rules! impl_vec_getter {
    (
        $type_name:path,
        $fn_name:ident,
        $ret_ty:ty,
        $helper:ident,
        |$obj:ident| $body:expr
    ) => {
        #[::safer_ffi::ffi_export]
        pub extern "C" fn $fn_name(
            h: ::std::option::Option<&$crate::core::ExprHandle>,
            out_err: ::std::option::Option<
                ::safer_ffi::prelude::Out<'_, $crate::core::TinnedErrorBox>,
            >,
        ) -> $ret_ty {
            $helper::<$type_name>(h, out_err, stringify!($fn_name), |$obj| $body)
        }
    };
}

macro_rules! impl_pert_multichain_getter {
    ($type_name:path, $fn_name:ident, |$obj:ident| $body:expr) => {
        #[::safer_ffi::ffi_export]
        pub extern "C" fn $fn_name(
            h: ::std::option::Option<&$crate::core::ExprHandle>,
            out_err: ::std::option::Option<
                ::safer_ffi::prelude::Out<'_, $crate::core::TinnedErrorBox>,
            >,
        ) -> ::std::option::Option<$crate::perturbations::PertMultichainBox> {
            $crate::c_support::ffi_map_expr_as::<$type_name, _>(
                h,
                out_err,
                stringify!($fn_name),
                |$obj| {
                    let chain = ::std::sync::Arc::new($body);
                    Ok($crate::perturbations::PertMultichainBox::new(
                        $crate::perturbations::PertMultichainHandle::new(chain),
                    ))
                },
            )
        }
    };
}
