// Compile-time test to make sure the struct implements Send + Sync
#[allow(unused_macros)]
macro_rules! test_struct_safety {
    ($type_name:ty) => {
        #[test]
        fn test_struct_send_sync() {
            fn assert_send_sync<T: Send + Sync>() {}
            assert_send_sync::<$type_name>();
        }
    };
}

// Test thread safety of interning
#[allow(unused_macros)]
macro_rules! test_thread_interning {
    ($make_expr:expr) => {
        #[test]
        fn test_thread_interning() {
            use std::thread;

            let mut handles = vec![];

            for _ in 0..10 {
                handles.push(thread::spawn(|| $make_expr));
            }

            let results: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();

            for i in 1..results.len() {
                assert!(Arc::ptr_eq(&results[0], &results[i]));
            }
        }
    };
}

// Test LagMultiplier, NonElecFunction, OneElecOperator and WfnParameter
#[allow(unused_macros)]
macro_rules! test_nullary_oper {
    ($type_name:ident, $oper_name:ident, $make_expr:ident, $has_deps:tt, $is_scalar:tt) => {
        test_struct_safety!($type_name);

        test_thread_interning!($make_expr($oper_name));

        #[test]
        fn test_impl_expr() {
            let deriv = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
            test_nullary_oper!(@test_nullary_expr
                $type_name,
                $oper_name,
                deriv,
                $has_deps,
                $is_scalar,
            );
        }

        #[test]
        fn test_serialization() {
            let op = $make_expr("");
            let json = serde_json::to_string(&op).unwrap();
            let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
            assert_eq!(&op, &deserialized);
        }

        #[test]
        fn test_utils() {
            let op1 = $make_expr($oper_name);
            let op2 = $make_expr($oper_name);
            let op3 = $make_expr("");

            assert!(Arc::ptr_eq(&op1, &op2));
            assert!(!Arc::ptr_eq(&op1, &op3));

            assert!(is_expr_type::<$type_name>(&op1));
            assert!(!is_zero_expr(&op1));
            assert!(!is_one_expr(&op1));
        }
    };

    (@test_nullary_expr
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

        assert!(is_zero_expr(&op0));

        let deps = make_super_multichain(&$deriv, 1u32);
        let op1 = $type_name::builder($oper_name)
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();

        let op = downcast_from_arc::<$type_name>(&op1).unwrap();
        assert_eq!(
            op,
            &$type_name {
                name: $oper_name.into(),
                dependencies: deps.clone(),
                derivative: $deriv.clone()
            }
        );

        assert_eq!(op.name(), $oper_name);
        assert_eq!(op.dependencies(), &deps);
        assert_eq!(op.derivative(), &$deriv);

        let op2 = op.builder_from($deriv.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        assert_eq!(
            op1.hash_key(),
            format!(
                "{}({}; [{}]; [{}])",
                stringify!($type_name),
                $oper_name,
                deps.hash_key(),
                $deriv.hash_key()
            )
        );
        assert_eq!(op1.is_scalar(), $is_scalar);
        assert_eq!(format!("{}", op1), format!("{}^({})", $oper_name, $deriv));

        let op3 = $type_name::builder($oper_name)
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();
        let op4 = $type_name::builder(random_alphanumeric($oper_name.len() as u32 + 1))
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

    (@test_nullary_expr
        $type_name:ident,
        $oper_name:ident,
        $deriv:ident,
        false,
        $is_scalar:tt,
    ) => {
        let op1 = $type_name::builder($oper_name).derivative($deriv.clone()).build().unwrap();

        let op = downcast_from_arc::<$type_name>(&op1).unwrap();
        assert_eq!(
            op,
            &$type_name {
                name: $oper_name.into(),
                derivative: $deriv.clone()
            }
        );

        assert_eq!(op.name(), $oper_name);
        assert_eq!(op.derivative(), &$deriv);

        let op2 = op.builder_from($deriv.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        assert_eq!(
            op1.hash_key(),
            format!("{}({}; [{}])", stringify!($type_name), $oper_name, $deriv.hash_key())
        );
        assert_eq!(op1.is_scalar(), $is_scalar);
        assert_eq!(format!("{}", op1), format!("{}^({})", $oper_name, $deriv));

        let op3 = $type_name::builder($oper_name).derivative($deriv.clone()).build().unwrap();
        let op4 = $type_name::builder(random_alphanumeric($oper_name.len() as u32 + 1))
            .derivative($deriv.clone())
            .build()
            .unwrap();
        let op5 = $type_name::builder($oper_name).build().unwrap();

        assert_eq!(&op1, &op3);
        assert_ne!(&op1, &op4);
        assert_ne!(&op1, &op5);
    };
}

// Test ExchCorrEnergy and ExchCorrPotential
#[allow(unused_macros)]
macro_rules! test_exch_corr {
    (
        $type_name:ident,
        $oper_name:ident,
        $make_expr:ident,
        $grid_expr_name:ident,
        $build_grid_expr:ident,
        $is_scalar:literal,
    ) => {
        test_struct_safety!($type_name);

        test_thread_interning!({
            $make_expr(
                $oper_name,
                Some(make_non_elec_function("weight")),
                Some(make_wfn_parameter("density")),
                Some(make_one_elec_operator("Omega")),
            )
        });

        #[test]
        fn test_impl_expr() {
            let weight = make_non_elec_function("");
            let density = make_wfn_parameter("");
            let overlap = make_one_elec_operator("");
            let op1 = $make_expr(
                $oper_name,
                Some(weight.clone()),
                Some(density.clone()),
                Some(overlap.clone()),
            );

            let op = downcast_from_arc::<$type_name>(&op1).unwrap();
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
                    derivative: PertMultichain::new(),
                }
            );

            assert_eq!(op.name(), $oper_name);
            assert_eq!(op.grid_weight(), &weight);
            assert_eq!(op.density_matrix(), &density);
            assert_eq!(op.overlap_distribution(), &overlap);
            assert_eq!(op.$grid_expr_name(), &$grid_expr_name);

            assert_eq!(
                op1.hash_key(),
                format!(
                    "{}({}; {}; {}; {}; [{}]; {})",
                    stringify!($type_name),
                    $oper_name,
                    weight.hash_key(),
                    density.hash_key(),
                    overlap.hash_key(),
                    PertMultichain::new().hash_key(),
                    $grid_expr_name.hash_key(),
                )
            );
            assert_eq!(op1.is_scalar(), $is_scalar);
            assert_eq!(format!("{}", op1), format!("{}[{}]", $oper_name, $grid_expr_name));

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

        #[test]
        fn test_serialization() {
            let op = $make_expr("", None, None, None);
            let json = serde_json::to_string(&op).unwrap();
            let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
            assert_eq!(&op, &deserialized);
        }

        #[test]
        fn test_utils() {
            let weight = make_non_elec_function("");
            let density = make_wfn_parameter("");
            let overlap = make_one_elec_operator("");
            let op1 = $make_expr(
                $oper_name,
                Some(weight.clone()),
                Some(density.clone()),
                Some(overlap.clone()),
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

            assert!(Arc::ptr_eq(&op1, &op2));
            assert!(!Arc::ptr_eq(&op1, &op3));
            assert!(!Arc::ptr_eq(&op1, &op4));
            assert!(!Arc::ptr_eq(&op1, &op5));
            assert!(!Arc::ptr_eq(&op1, &op6));

            assert!(is_expr_type::<$type_name>(&op1));
            assert!(!is_zero_expr(&op1));
            assert!(!is_one_expr(&op1));
        }
    };
}

// Test interning, serialization and utils for Trace, Transpose and HermitianTranspose
#[allow(unused_macros)]
macro_rules! test_unary_oper_properties {
    ($type_name:ident) => {
        test_struct_safety!($type_name);

        test_thread_interning!($type_name::new(make_two_elec_operator(
            "op(2el)",
            Some(make_wfn_parameter("wfn"))
        ))
        .unwrap());

        #[test]
        fn test_serialization() {
            let op = $type_name::new(make_two_elec_operator("", None)).unwrap();
            let json = serde_json::to_string(&op).unwrap();
            let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
            assert_eq!(&op, &deserialized);
        }

        #[test]
        fn test_utils() {
            let density = make_wfn_parameter("");
            let arg_2el = make_two_elec_operator("op(2el)", Some(density.clone()));
            let op1 = $type_name::new(arg_2el.clone()).unwrap();
            let op2 = $type_name::new(arg_2el).unwrap();
            let op3 = $type_name::new(make_two_elec_operator("", Some(density))).unwrap();
            let op4 = $type_name::new(make_two_elec_operator("op(2el)", None)).unwrap();

            assert!(Arc::ptr_eq(&op1, &op2));
            assert!(!Arc::ptr_eq(&op1, &op3));
            assert!(!Arc::ptr_eq(&op1, &op4));

            assert!(is_expr_type::<$type_name>(&op1));
            assert!(!is_zero_expr(&op1));
            assert!(!is_one_expr(&op1));
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
            let op0 = $type_name::new(ZeroOperator::new()).unwrap();
            assert!(is_zero_expr(&op0));

            let arg_2el = make_two_elec_operator("", None);
            let op1 = $type_name::new(arg_2el.clone()).unwrap();

            let op = downcast_from_arc::<$type_name>(&op1).unwrap();
            assert_eq!(
                op,
                &$type_name {
                    argument: arg_2el.clone()
                }
            );
            assert_eq!(op.argument(), &arg_2el);

            let op2 = $type_name::new(arg_2el.clone()).unwrap();
            assert!(Arc::ptr_eq(&op1, &op2));
            assert_eq!(&op1, &op2);

            assert_eq!(
                op1.hash_key(),
                format!("{}({})", stringify!($type_name), arg_2el.hash_key())
            );
            assert!(!op1.is_scalar());
            assert_eq!(format!("{}", op1), format!($display_fmt, arg = arg_2el));

            let mut argument = Conjugate::new(arg_2el.clone()).unwrap();
            let op3 = $type_name::new(argument).unwrap();
            assert_ne!(&op1, &op3);
            assert_eq!(&op3, &$conj_type::new(arg_2el.clone()).unwrap());

            let op4 = $type_name::new(op2).unwrap();
            assert_ne!(&op1, &op4);
            assert_eq!(&op4, &arg_2el);

            argument = $conj_type::new(arg_2el.clone()).unwrap();
            let op5 = $type_name::new(argument).unwrap();
            assert_ne!(&op1, &op5);
            assert_eq!(&op5, &Conjugate::new(arg_2el.clone()).unwrap());

            let coef = make_symbol(4u32);
            let arg_wfn = make_wfn_parameter("");
            argument = MatrixMul::new(vec![
                coef.clone(),
                arg_2el.clone(),
                arg_wfn.clone()
            ]).unwrap();
            let op6 = $type_name::new(argument).unwrap();
            assert_ne!(&op1, &op6);

            let matmul = downcast_from_arc::<MatrixMul>(&op6).unwrap();
            test_transpose!(@assert_matmul_coef matmul, coef, $has_conj);
            assert_eq!(
                matmul.factors(),
                vec![$type_name::new(MatrixMul::new(vec![arg_2el, arg_wfn]).unwrap()).unwrap()]
            );
        }
    };

    (@assert_matmul_coef $matmul:ident, $coef:ident, true) => {
        assert_eq!($matmul.coefficient(), &Conjugate::new($coef).unwrap())
    };

    (@assert_matmul_coef $matmul:ident, $coef:ident, false) => {
        assert_eq!($matmul.coefficient(), &$coef)
    };
}
