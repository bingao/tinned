use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::Number;

/// A scalar symbolic constant that becomes 0 after differentiation.
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Symbol {
    name: String,
}

impl Symbol {
    #[inline]
    pub fn new(name: impl Into<String>) -> Arc<dyn Expr> {
        crate::internal::intern_expr(Arc::new(Self {
            name: name.into(),
        }))
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }
}

impl ExprInternal for Symbol {
    impl_expr_internal_methods!(Symbol, false);

    #[inline]
    fn hash_key(&self) -> String {
        format!("Symbol({})", self.name)
    }

    #[inline]
    fn replace_one_in_children(
        &self,
        _expr: &Arc<dyn Expr>,
        _replacement: Arc<dyn Expr>,
        _include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(self.clone_expr())
    }

    #[inline]
    fn replace_all_in_children(
        &self,
        _map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        _include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(self.clone_expr())
    }
}

#[typetag::serde]
impl Expr for Symbol {
    impl_expr_common_methods!(true);

    #[inline]
    fn differentiate(
        &self,
        _s: Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(Number::zero())
    }

    #[inline]
    fn remove_one(&self, s: &Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_one_self(s, false) {
            Ok(Number::zero())
        } else {
            Ok(self.clone_expr())
        }
    }

    #[inline]
    fn remove_all(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_any_self(set, false) {
            Ok(Number::zero())
        } else {
            Ok(self.clone_expr())
        }
    }

    #[inline]
    fn retain_one(
        &self,
        s: &Arc<dyn Expr>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_one_self(s, include_derivatives) {
            Ok(self.clone_expr())
        } else {
            Ok(Number::zero())
        }
    }

    #[inline]
    fn retain_any(
        &self,
        set: &HashSet<Arc<dyn Expr>>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_any_self(set, include_derivatives) {
            Ok(self.clone_expr())
        } else {
            Ok(Number::zero())
        }
    }
}

impl std::fmt::Display for Symbol {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.name)
    }
}

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use rand::prelude::IndexedRandom;
    use rand::rng;

    #[inline]
    pub fn random_alphanumeric(len: u32) -> String {
        if len == 0 {
            "alpha".to_string()
        } else {
            let charset: &[u8] = b"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
            let mut rng = rng();

            (0..len)
                .map(|_| {
                    let c = charset.choose(&mut rng).unwrap();
                    *c as char
                })
                .collect()
        }
    }

    #[inline]
    pub fn make_symbol(len_name: u32) -> Arc<dyn Expr> {
        Symbol::new(random_alphanumeric(len_name))
    }
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;
    use crate::expressions::Add;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::{
        downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr, negate_expr, subtract_exprs,
    };

    test_struct_safety!(Symbol);

    test_thread_interning!(make_symbol(0u32));

    // Implementation for Expr
    #[test]
    fn test_impl_expr() {
        let s1 = make_symbol(0u32);
        let name = random_alphanumeric(0u32);

        let s = downcast_from_arc::<Symbol>(&s1).unwrap();
        assert_eq!(
            s,
            &Symbol {
                name: name.clone()
            }
        );
        assert_eq!(s.name(), name);

        assert_eq!(s1.hash_key(), format!("Symbol({name})"));
        assert!(s1.is_scalar());
        assert_eq!(format!("{}", s1), name);

        let s2 = make_symbol(0u32);
        let s3 = make_symbol(3u32);
        assert_eq!(&s1, &s2);
        assert_ne!(&s1, &s3);
    }

    #[test]
    fn test_differentiation() {
        let s = make_symbol(10u32);
        let p = make_perturbation_symbol(4u32, 4u32);
        assert_eq!(&s.differentiate(p).unwrap(), &Number::zero());
    }

    // Test serialization and deserialization via `serde_json`
    #[test]
    fn test_serialization() {
        let s = make_symbol(10u32);
        let json = serde_json::to_string(&s).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&s, &deserialized);
    }

    // Test utils
    #[test]
    fn test_utils() {
        let s1 = make_symbol(0u32);

        assert!(is_expr_type::<Symbol>(&s1));
        assert!(!is_zero_expr(&s1, None));
        assert!(!is_one_expr(&s1, None));

        let s2 = make_symbol(0u32);
        let s3 = make_symbol(10u32);

        assert!(Arc::ptr_eq(&s1, &s2));
        assert!(!Arc::ptr_eq(&s1, &s3));

        assert!(is_zero_expr(
            &Add::new(vec![s1.clone(), negate_expr(s1.clone()).unwrap()]).unwrap(),
            None
        ));
        assert!(is_zero_expr(&subtract_exprs(s1.clone(), s2.clone()).unwrap(), None));
    }
}
