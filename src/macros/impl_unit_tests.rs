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

        fn make_chain() -> PertMultichain {
            let mut map = BTreeMap::new();
            let p1 = Perturbation::new("p1", Symbol::new("omega"));
            let p2 = Perturbation::new("p2", Number::from_f64(3.14));

            map.insert(p1.clone(), 2);
            map.insert(p2.clone(), 4);

            let chain = PertMultichain::from_map(map);
            chain
        }

        test_thread_interning!(
            test_nullary_oper!(@make_nullary_expr $type_name, "op", $has_deps)
        );

        #[test]
        fn test_struct() {
            let deriv = make_chain();
            test_nullary_oper!(@test_nullary_struct $type_name, "op", deriv, $has_deps, $is_scalar);
        }

        #[test]
        fn test_impl_expr() {
            let deriv = make_chain();
            test_nullary_oper!(@test_nullary_expr $type_name, "op", deriv, $has_deps, $is_scalar);
        }

        #[test]
        fn test_serialization() {
            let op = test_nullary_oper!(@make_nullary_expr $type_name, "op", $has_deps);

            let json = serde_json::to_string(&op).unwrap();
            let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
            assert!(op == deserialized);
        }

        #[test]
        fn test_utils() {
            let deriv = make_chain();
            test_nullary_oper!(@test_nullary_utils $type_name, "op", deriv, $has_deps, $is_scalar);
        }
    };

    (@make_nullary_expr $type_name:ident, $oper_name:literal, true) => {{
        let deps = make_chain();
        let deriv = make_chain();
        $type_name::builder($oper_name).dependencies(deps).derivative(deriv).build().unwrap()
    }};

    (@make_nullary_expr $type_name:ident, $oper_name:literal, false) => {{
        let deriv = make_chain();
        $type_name::builder($oper_name).derivative(deriv).build().unwrap()
    }};

    (@test_nullary_struct
        $type_name:ident,
        $oper_name:literal,
        $deriv:ident,
        true,
        $is_scalar:tt
    ) => {
        let deps = make_chain();
        let op = $type_name {
            name: $oper_name.into(),
            dependencies: deps.clone(),
            derivative: $deriv.clone(),
        };

        assert_eq!(op.name(), $oper_name);
        assert_eq!(op.dependencies(), &deps.clone());
        assert_eq!(op.derivative(), &$deriv.clone());
        assert_eq!(
            op.hash_key(),
            format!(
                "{}({}; [{}]; [{}])",
                stringify!($type_name),
                $oper_name,
                deps.hash_key(),
                $deriv.hash_key()
            )
        );
        assert_eq!(op.is_scalar(), $is_scalar);

        assert_eq!(
            op,
            $type_name {
                name: $oper_name.into(),
                dependencies: deps.clone(),
                derivative: $deriv.clone(),
            }
        );
        assert_ne!(
            op,
            $type_name {
                name: "tinned".into(),
                dependencies: deps.clone(),
                derivative: $deriv.clone(),
            }
        );
        assert_ne!(
            op,
            $type_name {
                name: $oper_name.into(),
                dependencies: PertMultichain::new(),
                derivative: $deriv.clone(),
            }
        );
        assert_ne!(
            op,
            $type_name {
                name: $oper_name.into(),
                dependencies: deps.clone(),
                derivative: PertMultichain::new(),
            }
        );

        assert_eq!(format!("{}", op), format!("{}^({})", $oper_name, $deriv));
    };

    (@test_nullary_struct
        $type_name:ident,
        $oper_name:literal,
        $deriv:ident,
        false,
        $is_scalar:tt
    ) => {
        let op = $type_name {
            name: $oper_name.into(),
            derivative: $deriv.clone(),
        };

        assert_eq!(op.name(), $oper_name);
        assert_eq!(op.derivative(), &$deriv.clone());
        assert_eq!(
            op.hash_key(),
            format!("{}({}; [{}])", stringify!($type_name), $oper_name, $deriv.hash_key())
        );
        assert_eq!(op.is_scalar(), $is_scalar);

        assert_eq!(
            op,
            $type_name {
                name: $oper_name.into(),
                derivative: $deriv.clone(),
            }
        );
        assert_ne!(
            op,
            $type_name {
                name: "tinned".into(),
                derivative: $deriv.clone(),
            }
        );
        assert_ne!(
            op,
            $type_name {
                name: $oper_name.into(),
                derivative: PertMultichain::new(),
            }
        );

        assert_eq!(format!("{}", op), format!("{}^({})", $oper_name, $deriv));
    };

    (@test_nullary_expr
        $type_name:ident,
        $oper_name:literal,
        $deriv:ident,
        true,
        $is_scalar:tt
    ) => {
        let deps = make_chain();
        let op = $type_name::builder($oper_name)
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();

        assert_eq!(
            op.hash_key(),
            format!("{}({}; [{}]; [{}])", stringify!($type_name), $oper_name, deps, $deriv)
        );
        assert_eq!(op.is_scalar(), $is_scalar);
        assert_eq!(format!("{}", op), format!("{}^({})", $oper_name, $deriv));

        let op1 = $type_name::builder($oper_name)
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();
        let op2 = $type_name::builder("tinned")
            .dependencies(deps.clone())
            .derivative($deriv.clone())
            .build()
            .unwrap();
        let op3 = $type_name::builder($oper_name).dependencies(deps).build().unwrap();
        let op4 = $type_name::builder($oper_name).derivative($deriv).build().unwrap();

        assert!(op == op1);
        assert!(op != op2);
        assert!(op != op3);
        assert!(op != op4);
    };

    (@test_nullary_expr
        $type_name:ident,
        $oper_name:literal,
        $deriv:ident,
        false,
        $is_scalar:tt
    ) => {
        let op = $type_name::builder($oper_name).derivative($deriv.clone()).build().unwrap();

        assert_eq!(
            op.hash_key(),
            format!("{}({}; [{}])", stringify!($type_name), $oper_name, $deriv)
        );
        assert_eq!(op.is_scalar(), $is_scalar);
        assert_eq!(format!("{}", op), format!("{}^({})", $oper_name, $deriv));

        let op1 = $type_name::builder($oper_name).derivative($deriv.clone()).build().unwrap();
        let op2 = $type_name::builder("tinned").derivative($deriv.clone()).build().unwrap();
        let op3 = $type_name::builder($oper_name).build().unwrap();

        assert!(op == op1);
        assert!(op != op2);
        assert!(op != op3);
    };

    (@test_nullary_utils
        $type_name:ident,
        $oper_name:literal,
        $deriv:ident,
        true,
        $is_scalar:tt
    ) => {
        let deps = make_chain();
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

        let op = downcast_from_arc::<$type_name>(&op1).unwrap();
        assert_eq!(
            op,
            &$type_name {
                name: $oper_name.into(),
                dependencies: deps.clone(),
                derivative: $deriv.clone()
            }
        );

        assert!(is_expr_type::<$type_name>(&op1));
        assert!(!is_zero_expr(&op1));
        assert!(!is_one_expr(&op1));

        let op4 = op.builder_from($deriv).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op4));
        assert!(op1 == op4);
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

        let op = downcast_from_arc::<$type_name>(&op1).unwrap();
        assert_eq!(
            op,
            &$type_name {
                name: $oper_name.into(),
                derivative: $deriv.clone()
            }
        );

        assert!(is_expr_type::<$type_name>(&op1));
        assert!(!is_zero_expr(&op1));
        assert!(!is_one_expr(&op1));

        let op4 = op.builder_from($deriv).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op4));
        assert!(op1 == op4);
    };
}
