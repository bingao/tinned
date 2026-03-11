// Test ExchCorrEnergy and ExchCorrPotential
#[allow(unused_macros)]
macro_rules! test_exch_corr {
    (
        $type_name:ident,
        $oper_name:ident,
        $make_expr:ident,
        $grid_expr_name:ident,
        $build_grid_expr:ident,
        $is_scalar:literal
    ) => {
        test_struct_safety!($type_name);

        test_thread_interning!({
            $make_expr(
                $oper_name,
                Some($crate::expressions::non_elec_function::test_utils::make_non_elec_function(
                    "weight",
                )),
                Some($crate::expressions::wfn_parameter::test_utils::make_wfn_parameter("density")),
                Some($crate::expressions::one_elec_operator::test_utils::make_one_elec_operator(
                    "Omega",
                )),
            )
        });

        #[test]
        fn test_impl_expr() {
            let weight =
                $crate::expressions::non_elec_function::test_utils::make_non_elec_function("");
            let density = $crate::expressions::wfn_parameter::test_utils::make_wfn_parameter("");
            let overlap =
                $crate::expressions::one_elec_operator::test_utils::make_one_elec_operator("");
            let op1 = $make_expr(
                $oper_name,
                Some(weight.clone()),
                Some(density.clone()),
                Some(overlap.clone()),
            );

            let op = $crate::public::downcast_from_arc::<$type_name>(&op1).unwrap();
            let $grid_expr_name =
                $build_grid_expr(weight.clone(), density.clone(), overlap.clone()).unwrap();
            assert_eq!(
                op,
                &$type_name {
                    name: $oper_name.into(),
                    grid_weight: weight.clone(),
                    density_matrix: density.clone(),
                    overlap_distribution: overlap.clone(),
                    $grid_expr_name: $grid_expr_name.clone(),
                    derivative: $crate::perturbations::PertMultichain::new(),
                }
            );

            assert_eq!(op.name(), $oper_name);
            assert_eq!(op.grid_weight(), &weight);
            assert_eq!(op.density_matrix(), &density);
            assert_eq!(op.overlap_distribution(), &overlap);
            assert_eq!(op.$grid_expr_name(), &$grid_expr_name);

            assert_eq!(
                op1.hash_key(),
                ::std::format!(
                    "{}({}; {}; {}; {}; [{}]; {})",
                    stringify!($type_name),
                    $oper_name,
                    weight.hash_key(),
                    density.hash_key(),
                    overlap.hash_key(),
                    $crate::perturbations::PertMultichain::new().hash_key(),
                    $grid_expr_name.hash_key(),
                )
            );
            assert_eq!(op1.is_scalar(), $is_scalar);
            assert_eq!(
                ::std::format!("{}", op1),
                ::std::format!("{}[{}]", $oper_name, $grid_expr_name)
            );

            let op2 = $make_expr(
                $oper_name,
                Some(weight.clone()),
                Some(density.clone()),
                Some(overlap.clone()),
            );
            let op3 =
                $make_expr("", Some(weight.clone()), Some(density.clone()), Some(overlap.clone()));
            let op4 = $make_expr($oper_name, None, Some(density.clone()), Some(overlap.clone()));
            let op5 = $make_expr($oper_name, Some(weight.clone()), None, Some(overlap));
            let op6 = $make_expr($oper_name, Some(weight), Some(density), None);

            assert_eq!(&op1, &op2);
            assert_ne!(&op1, &op3);
            assert_ne!(&op1, &op4);
            assert_ne!(&op1, &op5);
            assert_ne!(&op1, &op6);
        }

        // This unit test will not check the correctness of differentiation on
        // XC energy or potenital at grid points, i.e. `$grid_expr_name`. That
        // will be checked in some integration tests.
        #[test]
        fn test_differentiation() {
            let weight =
                $crate::expressions::non_elec_function::test_utils::make_non_elec_function("");
            let density = $crate::expressions::wfn_parameter::test_utils::make_wfn_parameter("");
            let overlap =
                $crate::expressions::one_elec_operator::test_utils::make_one_elec_operator("");
            let op = $make_expr(
                $oper_name,
                Some(weight.clone()),
                Some(density.clone()),
                Some(overlap.clone()),
            );

            let mut p = make_perturbation_symbol(4u32, 4u32);
            let mut diff_op = op.differentiate(&p).unwrap();
            let mut deriv = $crate::perturbations::PertMultichain::new();
            deriv.insert(&p);

            let mut diff_cast = $crate::public::downcast_from_arc::<$type_name>(&diff_op).unwrap();

            assert_eq!(diff_cast.derivative(), &deriv);

            let $grid_expr_name =
                $build_grid_expr(weight.clone(), density.clone(), overlap.clone()).unwrap();
            let mut diff_grid_expr = $grid_expr_name.differentiate(&p).unwrap();

            assert_eq!(diff_cast.$grid_expr_name(), &diff_grid_expr);

            // The second order differentiation
            p = $crate::perturbations::perturbation::test_utils::make_perturbation_symbol(
                4u32, 4u32,
            );
            diff_op = diff_op.differentiate(&p).unwrap();
            deriv.insert(&p);
            diff_cast = $crate::public::downcast_from_arc::<$type_name>(&diff_op).unwrap();

            assert_eq!(diff_cast.derivative(), &deriv);

            diff_grid_expr = diff_grid_expr.differentiate(&p).unwrap();

            assert_eq!(diff_cast.$grid_expr_name(), &diff_grid_expr);
        }

        #[test]
        fn test_serialization() {
            let op = $make_expr("", None, None, None);
            let json = ::serde_json::to_string(&op).unwrap();
            let deserialized: expr_arc_ty!() = ::serde_json::from_str(&json).unwrap();
            assert_eq!(&op, &deserialized);
        }

        #[test]
        fn test_utils() {
            let weight =
                $crate::expressions::non_elec_function::test_utils::make_non_elec_function("");
            let density = $crate::expressions::wfn_parameter::test_utils::make_wfn_parameter("");
            let overlap =
                $crate::expressions::one_elec_operator::test_utils::make_one_elec_operator("");
            let op1 = $make_expr(
                $oper_name,
                Some(weight.clone()),
                Some(density.clone()),
                Some(overlap.clone()),
            );

            assert!($crate::public::is_expr_type::<$type_name>(&op1));
            assert!(!$crate::public::is_zero_expr(&op1, None));
            assert!(!$crate::public::is_one_expr(&op1, None));

            let op2 = $make_expr(
                $oper_name,
                Some(weight.clone()),
                Some(density.clone()),
                Some(overlap.clone()),
            );
            let op3 =
                $make_expr("", Some(weight.clone()), Some(density.clone()), Some(overlap.clone()));
            let op4 = $make_expr($oper_name, None, Some(density.clone()), Some(overlap.clone()));
            let op5 = $make_expr($oper_name, Some(weight.clone()), None, Some(overlap));
            let op6 = $make_expr($oper_name, Some(weight), Some(density), None);

            assert!(::std::sync::Arc::ptr_eq(&op1, &op2));
            assert!(!::std::sync::Arc::ptr_eq(&op1, &op3));
            assert!(!::std::sync::Arc::ptr_eq(&op1, &op4));
            assert!(!::std::sync::Arc::ptr_eq(&op1, &op5));
            assert!(!::std::sync::Arc::ptr_eq(&op1, &op6));
        }
    };
}
