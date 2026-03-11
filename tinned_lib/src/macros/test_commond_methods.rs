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
            let mut handles = ::std::vec::Vec::new();

            for _ in 0..10 {
                handles.push(::std::thread::spawn(|| $make_expr));
            }

            let results: ::std::vec::Vec<_> =
                handles.into_iter().map(|h| h.join().unwrap()).collect();

            for i in 1..results.len() {
                assert!(::std::sync::Arc::ptr_eq(&results[0], &results[i]));
            }
        }
    };
}
