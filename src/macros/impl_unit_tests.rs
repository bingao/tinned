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
    ($type_name:ident, $has_deps:tt, $is_scalar:tt) => {
        test_struct_safety!($type_name);

        test_thread_interning!(
            test_nullary_oper!(@make_nullary_expr $type_name, "op", 0u32, 0u32, 0u32, $has_deps)
        );

        #[test]
        fn test_impl_expr() {
            let deriv = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
            test_nullary_oper!(@test_nullary_expr $type_name, "op", deriv, $has_deps, $is_scalar);
        }

        #[test]
        fn test_serialization() {
            let op = test_nullary_oper!(@make_nullary_expr $type_name, "op", 2u32, 8u32, 10u32, $has_deps);

            let json = serde_json::to_string(&op).unwrap();
            let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
            assert!(op == deserialized);
        }

        #[test]
        fn test_utils() {
            let deriv = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
            test_nullary_oper!(@test_nullary_utils $type_name, "op", deriv, $has_deps, $is_scalar);
        }
    };

    (@make_nullary_expr
        $type_name:ident,
        $oper_name:literal,
        $len_name:literal,
        $val_range:literal,
        $order_range:literal,
        true
    ) => {{
        let deps = make_pert_multichain($len_name, $val_range, 1u32, $order_range);
        let deriv = make_pert_multichain($len_name, $val_range, 1u32, $order_range);
        $type_name::builder($oper_name).dependencies(deps).derivative(deriv).build().unwrap()
    }};

    (@make_nullary_expr
        $type_name:ident,
        $oper_name:literal,
        $len_name:literal,
        $val_range:literal,
        $order_range:literal,
        false
    ) => {{
        let deriv = make_pert_multichain($len_name, $val_range, 1u32, $order_range);
        $type_name::builder($oper_name).derivative(deriv).build().unwrap()
    }};

    (@test_nullary_expr
        $type_name:ident,
        $oper_name:literal,
        $deriv:ident,
        true,
        $is_scalar:tt
    ) => {
        let deps = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
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
        assert_eq!(op.dependencies(), &deps.clone());
        assert_eq!(op.derivative(), &$deriv.clone());

        let op2 = op.builder_from($deriv.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(op1 == op2);

        assert_eq!(
            op1.hash_key(),
            format!("{}({}; [{}]; [{}])", stringify!($type_name), $oper_name, deps, $deriv)
        );
        assert_eq!(op1.is_scalar(), $is_scalar);
        assert_eq!(format!("{}", op1), format!("{}^({})", $oper_name, $deriv));

        let op3 = $type_name::builder($oper_name)
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();
        let op4 = $type_name::builder("tinned")
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();
        let op5 = $type_name::builder($oper_name).dependencies(deps).build().unwrap();
        let op6 = $type_name::builder($oper_name).derivative($deriv).build().unwrap();

        assert!(op1 == op3);
        assert!(op1 != op4);
        assert!(op1 != op5);
        assert!(op1 != op6);
    };

    (@test_nullary_expr
        $type_name:ident,
        $oper_name:literal,
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
        assert_eq!(op.derivative(), &$deriv.clone());

        let op2 = op.builder_from($deriv.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(op1 == op2);

        assert_eq!(
            op1.hash_key(),
            format!("{}({}; [{}])", stringify!($type_name), $oper_name, $deriv)
        );
        assert_eq!(op1.is_scalar(), $is_scalar);
        assert_eq!(format!("{}", op1), format!("{}^({})", $oper_name, $deriv));

        let op3 = $type_name::builder($oper_name).derivative($deriv.clone()).build().unwrap();
        let op4 = $type_name::builder("tinned").derivative($deriv.clone()).build().unwrap();
        let op5 = $type_name::builder($oper_name).build().unwrap();

        assert!(op1 == op3);
        assert!(op1 != op4);
        assert!(op1 != op5);
    };

    (@test_nullary_utils
        $type_name:ident,
        $oper_name:literal,
        $deriv:ident,
        true,
        $is_scalar:tt
    ) => {
        let deps = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
        let op1 = $type_name::builder($oper_name)
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();
        let op2 = $type_name::builder($oper_name)
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();
        let op3 = $type_name::builder("tinned")
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));

        assert!(is_expr_type::<$type_name>(&op1));
        assert!(!is_zero_expr(&op1));
        assert!(!is_one_expr(&op1));
    };

    (@test_nullary_utils
        $type_name:ident,
        $oper_name:literal,
        $deriv:ident,
        false,
        $is_scalar:tt
    ) => {
        let op1 = $type_name::builder($oper_name).derivative($deriv.clone()).build().unwrap();
        let op2 = $type_name::builder($oper_name).derivative($deriv.clone()).build().unwrap();
        let op3 = $type_name::builder("tinned").derivative($deriv.clone()).build().unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));

        assert!(is_expr_type::<$type_name>(&op1));
        assert!(!is_zero_expr(&op1));
        assert!(!is_one_expr(&op1));
    };
}
