// Test LagMultiplier, NonElecFunction, OneElecOperator and WfnParameter
#[allow(unused_macros)]
macro_rules! test_nullary_expr {
    ($type_name:ident, $oper_name:ident, $make_expr:ident, $has_deps:tt, $is_scalar:tt) => {
        test_struct_safety!($type_name);

        test_thread_interning!($make_expr($oper_name));

        #[test]
        fn test_impl_expr() {
            let deriv = $crate::perturbations::pert_multichain::test_utils::make_pert_multichain(2u32, 8u32, 1u32, 10u32);
            test_nullary_expr!(
                @test_nullary_impl_expr
                $type_name,
                $oper_name,
                deriv,
                $has_deps,
                $is_scalar,
            );
        }

        test_nullary_expr!(@test_nullary_differentiation $type_name, $oper_name, $has_deps);

        #[test]
        fn test_serialization() {
            let op = $make_expr("");
            let json = ::serde_json::to_string(&op).unwrap();
            let deserialized: expr_arc_ty!() = ::serde_json::from_str(&json).unwrap();
            assert_eq!(&op, &deserialized);
        }

        #[test]
        fn test_utils() {
            let op1 = $make_expr($oper_name);

            assert!($crate::public::is_expr_type::<$type_name>(&op1));
            assert!(!$crate::public::is_zero_expr(&op1, None));
            assert!(!$crate::public::is_one_expr(&op1, None));

            let op2 = $make_expr($oper_name);
            let op3 = $make_expr("");

            assert!(::std::sync::Arc::ptr_eq(&op1, &op2));
            assert!(!::std::sync::Arc::ptr_eq(&op1, &op3));
        }
    };

    (@test_nullary_impl_expr
        $type_name:ident,
        $oper_name:ident,
        $deriv:ident,
        true,
        $is_scalar:tt,
    ) => {
        let op0 = $type_name::builder($oper_name)
            .derivative($deriv.clone())
            .build()
            .unwrap();

        assert!($crate::public::is_zero_expr(&op0, None));

        let deps = $crate::perturbations::pert_multichain::test_utils::make_super_multichain(&$deriv, 1u32);
        let op1 = $type_name::builder($oper_name)
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();

        let is_perturbing = true;

        let op = $crate::public::downcast_from_arc::<$type_name>(&op1).unwrap();
        assert_eq!(
            op,
            &$type_name {
                name: $oper_name.into(),
                dependencies: deps.clone(),
                derivative: $deriv.clone(),
                is_perturbing,
            }
        );

        assert_eq!(op.name(), $oper_name);
        assert_eq!(op.dependencies(), &deps);
        assert_eq!(op.derivative(), &$deriv);

        let op2 = op.with_derivative($deriv.clone()).build().unwrap();
        assert!(::std::sync::Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        assert_eq!(
            op1.hash_key(),
            ::std::format!(
                "{}({}; [{}]; [{}]; {})",
                stringify!($type_name),
                $oper_name,
                deps.hash_key(),
                $deriv.hash_key(),
                is_perturbing,
            )
        );
        assert_eq!(op1.is_scalar(), $is_scalar);
        assert_eq!(
            ::std::format!("{}", op1),
            ::std::format!("{}({})^({})", $oper_name, is_perturbing, $deriv)
        );

        let op3 = $type_name::builder($oper_name)
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();
        let op4 = $type_name::builder(
                $crate::expressions::symbol::test_utils::random_alphanumeric(
                    $oper_name.len() as u32 + 1
                )
            )
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();
        let op5 = $type_name::builder($oper_name).dependencies(deps).build().unwrap();
        let op6 = $type_name::builder($oper_name).derivative($deriv).build().unwrap();

        assert_eq!(&op1, &op3);
        assert_ne!(&op1, &op4);
        assert_ne!(&op1, &op5);
        assert_ne!(&op1, &op6);
    };

    (@test_nullary_impl_expr
        $type_name:ident,
        $oper_name:ident,
        $deriv:ident,
        false,
        $is_scalar:tt,
    ) => {
        let op1 = $type_name::builder($oper_name)
            .derivative($deriv.clone())
            .build()
            .unwrap();

        let op = $crate::public::downcast_from_arc::<$type_name>(&op1).unwrap();
        assert_eq!(
            op,
            &$type_name {
                name: $oper_name.into(),
                derivative: $deriv.clone()
            }
        );

        assert_eq!(op.name(), $oper_name);
        assert_eq!(op.derivative(), &$deriv);

        let op2 = op.with_derivative($deriv.clone()).build().unwrap();
        assert!(::std::sync::Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        assert_eq!(
            op1.hash_key(),
            ::std::format!(
                "{}({}; [{}])",
                stringify!($type_name),
                $oper_name,
                $deriv.hash_key()
            )
        );
        assert_eq!(op1.is_scalar(), $is_scalar);
        assert_eq!(
            ::std::format!("{}", op1),
            ::std::format!("{}^({})", $oper_name, $deriv)
        );

        let op3 = $type_name::builder($oper_name)
            .derivative($deriv.clone())
            .build()
            .unwrap();
        let op4 = $type_name::builder(
                $crate::expressions::symbol::test_utils::random_alphanumeric(
                    $oper_name.len() as u32 + 1
                )
            )
            .derivative($deriv.clone())
            .build()
            .unwrap();
        let op5 = $type_name::builder($oper_name).build().unwrap();

        assert_eq!(&op1, &op3);
        assert_ne!(&op1, &op4);
        assert_ne!(&op1, &op5);
    };

    (@test_nullary_differentiation $type_name:ident, $oper_name:ident, true) => {
        #[test]
        fn test_differentiation() {
            let len_pert_name: u32 = 2;
            let mut deriv = $crate::perturbations::pert_multichain::test_utils::make_pert_multichain(len_pert_name, 8u32, 1u32, 10u32);
            let deps = $crate::perturbations::pert_multichain::test_utils::make_super_multichain(&deriv, 1u32);
            let op = $type_name::builder($oper_name)
                .dependencies(deps.clone())
                .derivative(deriv.clone())
                .build()
                .unwrap();

            let p: ::std::sync::Arc<$crate::perturbations::Perturbation> =
                deps.keys().first().cloned().unwrap();
            let diff_op = op.differentiate(&p).unwrap();
            deriv.insert(&p);

            let diff_cast = $crate::public::downcast_from_arc::<$type_name>(&diff_op).unwrap();
            assert_eq!(diff_cast.derivative(), &deriv);

            assert!($crate::public::is_zero_expr(
                &diff_op.differentiate(&p).unwrap(),
                None
            ));

            assert!($crate::public::is_zero_expr(
                &op.differentiate(&$crate::perturbations::perturbation::test_utils::make_perturbation_symbol(len_pert_name + 1u32, 4u32))
                    .unwrap(),
                None
            ));
        }
    };

    (@test_nullary_differentiation $type_name:ident, $oper_name:ident, false) => {
        #[test]
        fn test_differentiation() {
            let len_pert_name: u32 = 2;
            let mut deriv = $crate::perturbations::pert_multichain::test_utils::make_pert_multichain(len_pert_name, 8u32, 1u32, 10u32);
            let op = $type_name::builder($oper_name)
                .derivative(deriv.clone())
                .build()
                .unwrap();

            let p = $crate::perturbations::perturbation::test_utils::make_perturbation_symbol(len_pert_name + 1u32, 4u32);
            let diff_op = op.differentiate(&p).unwrap();
            deriv.insert(&p);

            let diff_cast = $crate::public::downcast_from_arc::<$type_name>(&diff_op).unwrap();
            assert_eq!(diff_cast.derivative(), &deriv);
        }
    };
}
