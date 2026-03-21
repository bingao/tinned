macro_rules! impl_exch_corr_ffi {
    ($type_name:path, $stem:ident, $grid_expr_name:ident) => {
        ::paste::paste! {
            #[::safer_ffi::ffi_export]
            pub extern "C" fn [<tinned_ $stem _new>](
                name: ::std::option::Option<::safer_ffi::prelude::char_p::Ref<'_>>,
                grid_weight: ::std::option::Option<&$crate::core::ExprHandle>,
                density_matrix: ::std::option::Option<&$crate::core::ExprHandle>,
                overlap_distribution: ::std::option::Option<&$crate::core::ExprHandle>,
                out_err: ::std::option::Option<
                    ::safer_ffi::prelude::Out<'_, $crate::core::TinnedErrorBox>
                >,
            ) -> ::std::option::Option<$crate::core::ExprBox> {
                let Some(name) = $crate::c_support::tinned_string_from_cstr(name) else {
                    $crate::core::tinned_error_new(
                        out_err,
                        ::tinned::public::generic_error(
                            concat!(
                                "Null or invalid name passed to ",
                                stringify!([<tinned_ $stem _new>])
                            ),
                            None,
                        ),
                    );
                    return None;
                };

                let Some(grid_weight) = grid_weight else {
                    $crate::core::tinned_error_new(
                        out_err,
                        ::tinned::public::generic_error(
                            concat!(
                                "Null grid weight passed to ",
                                stringify!([<tinned_ $stem _new>])
                            ),
                            None,
                        ),
                    );
                    return None;
                };
                let grid_weight_arc = grid_weight.clone_arc();

                let Some(density_matrix) = density_matrix else {
                    $crate::core::tinned_error_new(
                        out_err,
                        ::tinned::public::generic_error(
                            concat!(
                                "Null density matrix passed to ",
                                stringify!([<tinned_ $stem _new>])
                            ),
                            None,
                        ),
                    );
                    return None;
                };
                let density_matrix_arc = density_matrix.clone_arc();

                let Some(overlap_distribution) = overlap_distribution else {
                    $crate::core::tinned_error_new(
                        out_err,
                        ::tinned::public::generic_error(
                            concat!(
                                "Null overlap distribution passed to ",
                                stringify!([<tinned_ $stem _new>])
                            ),
                            None,
                        ),
                    );
                    return None;
                };
                let overlap_distribution_arc = overlap_distribution.clone_arc();

                match $type_name::builder(
                    name,
                    grid_weight_arc,
                    density_matrix_arc,
                    overlap_distribution_arc,
                )
                .build() {
                    Ok(expr_arc) => Some(
                        $crate::core::ExprBox::new(
                            $crate::core::ExprHandle::new(expr_arc)
                        )
                    ),
                    Err(e) => {
                        $crate::core::tinned_error_new(out_err, e);
                        None
                    },
                }
            }

            impl_cstr_getter!(
                [<tinned_ $stem _name>] : $type_name => |xc| xc.name().to_string()
            );

            impl_expr_getters!(
                $type_name;
                [<tinned_ $stem _grid_weight>] => |xc| Ok(::std::sync::Arc::clone(xc.grid_weight())),
                [<tinned_ $stem _density_matrix>] => |xc| Ok(::std::sync::Arc::clone(xc.density_matrix())),
                [<tinned_ $stem _overlap_distribution>] => |xc| Ok(::std::sync::Arc::clone(xc.overlap_distribution())),
                [<tinned_ $stem _ $grid_expr_name>] => |xc| Ok(::std::sync::Arc::clone(xc.$grid_expr_name())),
            );

            impl_pert_multichain_getter!(
                $type_name, [<tinned_ $stem _derivative>], |xc| xc.derivative().clone()
            );
        }
    };
}
