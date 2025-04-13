use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::perturbations::{PertMultichain, Perturbation};

impl_nullary_oper_type!(NonElecFunction, NonElecFunctionBuilder, true, true);
impl_nullary_oper_traits!(NonElecFunction, true, true);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::{Number, Symbol};
    use crate::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_struct_safety!(NonElecFunction);

    test_thread_interning!({
        let f1: Arc<dyn Expr> = Number::Integer(1).into();
        let f2: Arc<dyn Expr> = Symbol::new("w");
        let mut chain = PertMultichain::new();
        chain.insert(&Perturbation::new("alpha", f1.clone()));
        chain.insert(&Perturbation::new("beta", f2.clone()));
        NonElecFunction::builder("hnuc").derivative(chain).build().unwrap()
    });

    #[test]
    fn test_struct() {
        let f1: Arc<dyn Expr> = Number::Integer(1).into();
        let f2: Arc<dyn Expr> = Symbol::new("w");
        let mut chain = PertMultichain::new();
        chain.insert(&Perturbation::new("alpha", f1.clone()));
        chain.insert(&Perturbation::new("beta", f2.clone()));

        let hnuc = NonElecFunction {
            name: "hnuc".into(),
            dependencies: PertMultichain::new(),
            derivative: chain.clone(),
        };

        assert_eq!(hnuc.name(), "hnuc");
        assert_eq!(hnuc.derivative(), &chain.clone());
        assert_eq!(hnuc.hash_key(), format!("NonElecFunction(hnuc; []; [{chain}])"));
        assert!(hnuc.is_scalar());

        assert_eq!(
            hnuc,
            NonElecFunction {
                name: "hnuc".into(),
                dependencies: PertMultichain::new(),
                derivative: chain.clone()
            }
        );
        assert_ne!(
            hnuc,
            NonElecFunction {
                name: "vnuc".into(),
                dependencies: PertMultichain::new(),
                derivative: chain.clone()
            }
        );
        assert_ne!(
            hnuc,
            NonElecFunction {
                name: "hnuc".into(),
                dependencies: PertMultichain::new(),
                derivative: PertMultichain::new()
            }
        );

        assert_eq!(format!("{}", hnuc), format!("hnuc^({chain})"));
    }

    #[test]
    fn test_impl_expr() {
        let f1: Arc<dyn Expr> = Number::Integer(1).into();
        let f2: Arc<dyn Expr> = Symbol::new("w");
        let mut chain = PertMultichain::new();
        chain.insert(&Perturbation::new("alpha", f1.clone()));
        chain.insert(&Perturbation::new("beta", f2.clone()));

        let hnuc = NonElecFunction::builder("hnuc").derivative(chain.clone()).build().unwrap();

        assert_eq!(hnuc.hash_key(), format!("NonElecFunction(hnuc; []; [{chain}])"));
        assert!(hnuc.is_scalar());
        assert_eq!(format!("{}", hnuc), format!("hnuc^({chain})"));

        let hnuc1 = NonElecFunction::builder("hnuc").derivative(chain.clone()).build().unwrap();
        let hnuc2 = NonElecFunction::builder("vnuc").derivative(chain.clone()).build().unwrap();

        assert!(hnuc == hnuc1);
        assert!(hnuc != hnuc2);
    }

    #[test]
    fn test_serialization() {
        let f1: Arc<dyn Expr> = Number::Integer(1).into();
        let f2: Arc<dyn Expr> = Symbol::new("w");
        let mut chain = PertMultichain::new();
        chain.insert(&Perturbation::new("alpha", f1.clone()));
        chain.insert(&Perturbation::new("beta", f2.clone()));

        let hnuc = NonElecFunction::builder("hnuc").derivative(chain).build().unwrap();

        let json = serde_json::to_string(&hnuc).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert!(hnuc == deserialized);
    }

    #[test]
    fn test_utils() {
        let f1: Arc<dyn Expr> = Number::Integer(1).into();
        let f2: Arc<dyn Expr> = Symbol::new("w");
        let mut chain = PertMultichain::new();
        chain.insert(&Perturbation::new("alpha", f1.clone()));
        chain.insert(&Perturbation::new("beta", f2.clone()));

        let hnuc1 = NonElecFunction::builder("hnuc").derivative(chain.clone()).build().unwrap();
        let hnuc2 = NonElecFunction::builder("hnuc").derivative(chain.clone()).build().unwrap();
        let hnuc3 = NonElecFunction::builder("vnuc").derivative(chain.clone()).build().unwrap();

        assert!(Arc::ptr_eq(&hnuc1, &hnuc2));
        assert!(!Arc::ptr_eq(&hnuc1, &hnuc3));

        let hnuc = downcast_from_arc::<NonElecFunction>(&hnuc1).unwrap();
        assert_eq!(
            hnuc,
            &NonElecFunction {
                name: "hnuc".into(),
                dependencies: PertMultichain::new(),
                derivative: chain.clone()
            }
        );

        assert!(is_expr_type::<NonElecFunction>(&hnuc1));
        assert!(!is_zero_expr(&hnuc1));
        assert!(!is_one_expr(&hnuc1));

        let hnuc4 = hnuc.builder_from(chain.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&hnuc1, &hnuc4));
        assert!(hnuc1 == hnuc4);
    }
}
