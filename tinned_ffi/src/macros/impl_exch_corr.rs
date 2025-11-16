macro_rules! impl_exch_corr_ffi {
    ($type_name:path, $stem:ident, $grid_expr_name:ident) => {
        paste! {
            #[ffi_export]
            pub extern "C" fn [<tinned_ $stem _new>](
                name: Option<char_p::Ref<'_>>,
                grid_weight: Option<&ExprHandle>,
                density_matrix: Option<&ExprHandle>,
                overlap_distribution: Option<&ExprHandle>,
                out_err: Option<Out<'_, TinnedErrorBox>>,
            ) -> Option<ExprBox> {
                let Some(name) = tinned_string_from_cstr(name) else {
                    tinned_error_new(
                        out_err,
                        generic_error(concat!("Null or invalid name passed to ", stringify!([<tinned_ $stem _new>])), None),
                    );
                    return None;
                };

                let Some(grid_weight) = grid_weight else {
                    tinned_error_new(
                        out_err,
                        generic_error(concat!("Null grid weight passed to ", stringify!([<tinned_ $stem _new>])), None),
                    );
                    return None;
                };
                let grid_weight_arc = grid_weight.clone_arc();

                let Some(density_matrix) = density_matrix else {
                    tinned_error_new(
                        out_err,
                        generic_error(concat!("Null density matrix passed to ", stringify!([<tinned_ $stem _new>])), None),
                    );
                    return None;
                };
                let density_matrix_arc = density_matrix.clone_arc();

                let Some(overlap_distribution) = overlap_distribution else {
                    tinned_error_new(
                        out_err,
                        generic_error(concat!("Null overlap distribution passed to ", stringify!([<tinned_ $stem _new>])), None),
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
                    Ok(expr_arc) => Some(ExprBox::new(ExprHandle::new(expr_arc))),
                    Err(e) => {
                        tinned_error_new(out_err, e);
                        None
                    },
                }
            }

            impl_cstr_getter!(
                [<tinned_ $stem _name>] : $type_name => |xc| xc.name().to_string()
            );

            impl_expr_getters!(
                $type_name;
                [<tinned_ $stem _grid_weight>] => |xc| Ok(Arc::clone(xc.grid_weight())),
                [<tinned_ $stem _density_matrix>] => |xc| Ok(Arc::clone(xc.density_matrix())),
                [<tinned_ $stem _overlap_distribution>] => |xc| Ok(Arc::clone(xc.overlap_distribution())),
                [<tinned_ $stem _ $grid_expr_name>] => |xc| Ok(Arc::clone(xc.$grid_expr_name())),
            );

            impl_pert_multichain_getter!(
                [<tinned_ $stem _derivative>] : $type_name => |xc| xc.derivative().clone()
            );
        }
    };
}
