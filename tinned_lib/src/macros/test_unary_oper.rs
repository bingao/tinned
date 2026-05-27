// Test common properties for `Trace`, `Transpose`, such as interning,
// differentiation, serialization and utils, etc.
#[allow(unused_macros)]
macro_rules! test_unary_oper_properties {
    ($type_name:ident, $build_expr:expr) => {
        test_struct_safety!($type_name);

        test_thread_interning!(
            ($build_expr)(
                $crate::expressions::ao_two_elec_matrix::test_utils::make_ao_two_elec_matrix(
                    "G^{AO}",
                    Some($crate::expressions::wfn_parameter::test_utils::make_wfn_parameter("wfn"))
                )
            )
            .unwrap()
        );

        #[test]
        fn test_differentiation() {
            let argument =
                $crate::expressions::ao_two_elec_matrix::test_utils::make_ao_two_elec_matrix(
                    "", None,
                );
            let op = ($build_expr)(argument.clone()).unwrap();
            let p = $crate::perturbations::perturbation::test_utils::make_perturbation_symbol(
                4u32, 4u32,
            );
            let diff_op = &op.differentiate(p.clone()).unwrap();
            let diff_cast = $crate::public::downcast_from_arc::<$type_name>(&diff_op).unwrap();

            assert_eq!(diff_cast.argument(), &argument.differentiate(p).unwrap());
        }

        #[test]
        fn test_serialization() {
            let op = ($build_expr)(
                $crate::expressions::ao_two_elec_matrix::test_utils::make_ao_two_elec_matrix(
                    "", None,
                ),
            )
            .unwrap();
            let json = ::serde_json::to_string(&op).unwrap();
            let deserialized: expr_arc_ty!() = ::serde_json::from_str(&json).unwrap();
            assert_eq!(&op, &deserialized);
        }

        #[test]
        fn test_utils() {
            let density = $crate::expressions::wfn_parameter::test_utils::make_wfn_parameter("");
            let arg_2el =
                $crate::expressions::ao_two_elec_matrix::test_utils::make_ao_two_elec_matrix(
                    "G^{AO}",
                    Some(density.clone()),
                );
            let op1 = ($build_expr)(arg_2el.clone()).unwrap();

            assert!($crate::public::is_expr_type::<$type_name>(&op1));
            assert!(!$crate::public::is_zero_expr(&op1, None));
            assert!(!$crate::public::is_one_expr(&op1, None));

            let op2 = ($build_expr)(arg_2el).unwrap();
            let op3 = ($build_expr)(
                $crate::expressions::ao_two_elec_matrix::test_utils::make_ao_two_elec_matrix(
                    "",
                    Some(density),
                ),
            )
            .unwrap();
            let op4 = ($build_expr)(
                $crate::expressions::ao_two_elec_matrix::test_utils::make_ao_two_elec_matrix(
                    "G^{AO}", None,
                ),
            )
            .unwrap();

            assert!(::std::sync::Arc::ptr_eq(&op1, &op2));
            assert!(!::std::sync::Arc::ptr_eq(&op1, &op3));
            assert!(!::std::sync::Arc::ptr_eq(&op1, &op4));
        }
    };
}
