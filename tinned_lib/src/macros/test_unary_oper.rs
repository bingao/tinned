// Test common properties for Trace, Transpose and HermitianTranspose, such as
// interning, differentiation, serialization and utils, etc.
#[allow(unused_macros)]
macro_rules! test_unary_oper_properties {
    ($type_name:ident) => {
        test_struct_safety!($type_name);

        test_thread_interning!(
            $type_name::new(
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
            let op = $type_name::new(argument.clone()).unwrap();
            let p = $crate::perturbations::perturbation::test_utils::make_perturbation_symbol(
                4u32, 4u32,
            );
            let diff_op = &op.differentiate(&p).unwrap();
            let diff_cast = $crate::public::downcast_from_arc::<$type_name>(&diff_op).unwrap();

            assert_eq!(diff_cast.argument(), &argument.differentiate(&p).unwrap());
        }

        #[test]
        fn test_serialization() {
            let op = $type_name::new(
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
            let op1 = $type_name::new(arg_2el.clone()).unwrap();

            assert!($crate::public::is_expr_type::<$type_name>(&op1));
            assert!(!$crate::public::is_zero_expr(&op1, None));
            assert!(!$crate::public::is_one_expr(&op1, None));

            let op2 = $type_name::new(arg_2el).unwrap();
            let op3 = $type_name::new(
                $crate::expressions::ao_two_elec_matrix::test_utils::make_ao_two_elec_matrix(
                    "",
                    Some(density),
                ),
            )
            .unwrap();
            let op4 = $type_name::new(
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

// Test Transpose and HermitianTranspose
#[allow(unused_macros)]
macro_rules! test_transpose {
    ($type_name:ident, $conj_type:ident, $has_conj:tt, $display_fmt:expr) => {
        test_unary_oper_properties!($type_name);

        #[test]
        fn test_impl_expr() {
            let op0 = $type_name::new($crate::expressions::ZeroOperator::new()).unwrap();
            assert!($crate::public::is_zero_expr(&op0, None));

            let arg_2el = $crate::expressions::ao_two_elec_matrix::test_utils::make_ao_two_elec_matrix("", None);
            let op1 = $type_name::new(arg_2el.clone()).unwrap();

            let op = $crate::public::downcast_from_arc::<$type_name>(&op1).unwrap();
            assert_eq!(
                op,
                &$type_name {
                    argument: arg_2el.clone()
                }
            );
            assert_eq!(op.argument(), &arg_2el);

            let op2 = $type_name::new(arg_2el.clone()).unwrap();
            assert!(::std::sync::Arc::ptr_eq(&op1, &op2));
            assert_eq!(&op1, &op2);

            assert_eq!(
                op1.hash_key(),
                ::std::format!("{}({})", stringify!($type_name), arg_2el.hash_key())
            );
            assert!(!op1.is_scalar());
            assert_eq!(
                ::std::format!("{}", op1),
                ::std::format!($display_fmt, arg = arg_2el)
            );

            let mut argument = $crate::expressions::Conjugate::new(arg_2el.clone()).unwrap();
            let op3 = $type_name::new(argument).unwrap();
            assert_ne!(&op1, &op3);
            assert_eq!(&op3, &$conj_type::new(arg_2el.clone()).unwrap());

            let op4 = $type_name::new(op2).unwrap();
            assert_ne!(&op1, &op4);
            assert_eq!(&op4, &arg_2el);

            argument = $conj_type::new(arg_2el.clone()).unwrap();
            let op5 = $type_name::new(argument).unwrap();
            assert_ne!(&op1, &op5);
            assert_eq!(
                &op5,
                &$crate::expressions::Conjugate::new(arg_2el.clone()).unwrap()
            );

            let coef = $crate::expressions::symbol::test_utils::make_symbol(4u32);
            let arg_wfn =
                $crate::expressions::wfn_parameter::test_utils::make_wfn_parameter("");
            argument = $crate::expressions::MatrixMul::new(::std::vec![
                coef.clone(),
                arg_2el.clone(),
                arg_wfn.clone()
            ])
            .unwrap();
            let op6 = $type_name::new(argument).unwrap();
            assert_ne!(&op1, &op6);

            let mat_mul =
                $crate::public::downcast_from_arc::<$crate::expressions::MatrixMul>(&op6).unwrap();
            test_transpose!(@assert_mat_mul_coef mat_mul, coef, $has_conj);
            assert_eq!(
                mat_mul.factors(),
                ::std::vec![
                    $type_name::new(
                        $crate::expressions::MatrixMul::new(::std::vec![arg_2el, arg_wfn])
                            .unwrap()
                    )
                    .unwrap()
                ]
            );
        }
    };

    (@assert_mat_mul_coef $mat_mul:ident, $coef:ident, true) => {
        assert_eq!(
            $mat_mul.coefficient(),
            &$crate::expressions::Conjugate::new($coef).unwrap()
        )
    };

    (@assert_mat_mul_coef $mat_mul:ident, $coef:ident, false) => {
        assert_eq!($mat_mul.coefficient(), &$coef)
    };
}
