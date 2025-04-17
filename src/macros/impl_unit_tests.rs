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
            test_nullary_oper!(@test_nullary_expr $type_name, $oper_name, deriv, $has_deps, $is_scalar);
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
        $is_scalar:tt
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
        $is_scalar:tt
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
