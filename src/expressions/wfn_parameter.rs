use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::perturbations::{PertMultichain, Perturbation};

impl_nullary_oper_type!(WfnParameter, WfnParameterBuilder, false, false);
impl_nullary_oper_traits!(WfnParameter, false, false);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::{Number, Symbol};
    use crate::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    test_struct_safety!(WfnParameter);

    test_thread_interning!({
        let f1: Arc<dyn Expr> = Number::Integer(1).into();
        let f2: Arc<dyn Expr> = Symbol::new("w");
        let mut chain = PertMultichain::new();
        chain.insert(&Perturbation::new("alpha", f1.clone()));
        chain.insert(&Perturbation::new("beta", f2.clone()));
        WfnParameter::builder("Wfn").derivative(chain).build().unwrap()
    });

    #[test]
    fn test_struct() {
        let f1: Arc<dyn Expr> = Number::Integer(1).into();
        let f2: Arc<dyn Expr> = Symbol::new("w");
        let mut chain = PertMultichain::new();
        chain.insert(&Perturbation::new("alpha", f1.clone()));
        chain.insert(&Perturbation::new("beta", f2.clone()));

        let wfn = WfnParameter {
            name: "Wfn".into(),
            derivative: chain.clone(),
        };

        assert_eq!(wfn.name(), "Wfn");
        assert_eq!(wfn.derivative(), &chain.clone());
        assert_eq!(wfn.hash_key(), format!("WfnParameter(Wfn; [{chain}])"));
        assert!(!wfn.is_scalar());

        assert_eq!(
            wfn,
            WfnParameter {
                name: "Wfn".into(),
                derivative: chain.clone()
            }
        );
        assert_ne!(
            wfn,
            WfnParameter {
                name: "Other".into(),
                derivative: chain.clone()
            }
        );
        assert_ne!(
            wfn,
            WfnParameter {
                name: "Wfn".into(),
                derivative: PertMultichain::new()
            }
        );

        assert_eq!(format!("{}", wfn), format!("Wfn^({chain})"));
    }

    #[test]
    fn test_impl_expr() {
        let f1: Arc<dyn Expr> = Number::Integer(1).into();
        let f2: Arc<dyn Expr> = Symbol::new("w");
        let mut chain = PertMultichain::new();
        chain.insert(&Perturbation::new("alpha", f1.clone()));
        chain.insert(&Perturbation::new("beta", f2.clone()));

        let wfn = WfnParameter::builder("Wfn").derivative(chain.clone()).build().unwrap();

        assert_eq!(wfn.hash_key(), format!("WfnParameter(Wfn; [{chain}])"));
        assert!(!wfn.is_scalar());
        assert_eq!(format!("{}", wfn), format!("Wfn^({chain})"));

        let wfn1 = WfnParameter::builder("Wfn").derivative(chain.clone()).build().unwrap();
        let wfn2 = WfnParameter::builder("Other").derivative(chain.clone()).build().unwrap();

        assert!(wfn == wfn1);
        assert!(wfn != wfn2);
    }

    #[test]
    fn test_serialization() {
        let f1: Arc<dyn Expr> = Number::Integer(1).into();
        let f2: Arc<dyn Expr> = Symbol::new("w");
        let mut chain = PertMultichain::new();
        chain.insert(&Perturbation::new("alpha", f1.clone()));
        chain.insert(&Perturbation::new("beta", f2.clone()));

        let wfn = WfnParameter::builder("Wfn").derivative(chain).build().unwrap();

        let json = serde_json::to_string(&wfn).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert!(wfn == deserialized);
    }

    #[test]
    fn test_utils() {
        let f1: Arc<dyn Expr> = Number::Integer(1).into();
        let f2: Arc<dyn Expr> = Symbol::new("w");
        let mut chain = PertMultichain::new();
        chain.insert(&Perturbation::new("alpha", f1.clone()));
        chain.insert(&Perturbation::new("beta", f2.clone()));

        let wfn1 = WfnParameter::builder("Wfn").derivative(chain.clone()).build().unwrap();
        let wfn2 = WfnParameter::builder("Wfn").derivative(chain.clone()).build().unwrap();
        let wfn3 = WfnParameter::builder("Other").derivative(chain.clone()).build().unwrap();

        assert!(Arc::ptr_eq(&wfn1, &wfn2));
        assert!(!Arc::ptr_eq(&wfn1, &wfn3));

        let wfn = downcast_from_arc::<WfnParameter>(&wfn1).unwrap();
        assert_eq!(
            wfn,
            &WfnParameter {
                name: "Wfn".into(),
                derivative: chain.clone()
            }
        );

        assert!(is_expr_type::<WfnParameter>(&wfn1));
        assert!(!is_zero_expr(&wfn1));
        assert!(!is_one_expr(&wfn1));

        let wfn4 = wfn.builder_from(chain.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&wfn1, &wfn4));
        assert!(wfn1 == wfn4);
    }
}
