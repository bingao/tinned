use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::ZeroOperator;

/// Excitation operator, which can not be differentiated. Implemented as a
/// non-scalar `Symbol`.
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ExcitationOperator {
    name: String,
}

impl ExcitationOperator {
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

impl ExprInternal for ExcitationOperator {
    impl_expr_internal_methods!(ExcitationOperator, false);

    #[inline]
    fn hash_key(&self) -> String {
        format!("ExcitationOperator({})", self.name)
    }

    #[inline]
    fn replace_expr_fields(
        &self,
        _map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        _exact_equality: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(self.clone_expr())
    }

    // `retain_expr_fields` may be called by `MatrixMul`
    #[inline]
    fn retain_expr_fields(
        &self,
        _expr: &Arc<dyn Expr>,
        _exact_equality: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(ZeroOperator::new())
    }
}

#[typetag::serde]
impl Expr for ExcitationOperator {
    impl_expr_common_methods!(false);

    #[inline]
    fn differentiate(
        &self,
        _s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(ZeroOperator::new())
    }

    #[inline]
    fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
            Ok(ZeroOperator::new())
        } else {
            Ok(self.clone_expr())
        }
    }

    #[inline]
    fn retain(
        &self,
        set: &HashSet<Arc<dyn Expr>>,
        _exact_equality: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if set.iter().all(|expr| self.eq_expr(expr.as_ref())) {
            Ok(self.clone_expr())
        } else {
            Ok(ZeroOperator::new())
        }
    }
}

impl std::fmt::Display for ExcitationOperator {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.name)
    }
}

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use crate::expressions::symbol::test_utils::random_alphanumeric;

    #[inline]
    pub fn make_excitation_operator(len_name: u32) -> Arc<dyn Expr> {
        ExcitationOperator::new(random_alphanumeric(len_name))
    }
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;
    use crate::expressions::MatrixAdd;
    use crate::expressions::symbol::test_utils::random_alphanumeric;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::{
        downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr, negate_expr, subtract_exprs,
    };

    test_struct_safety!(ExcitationOperator);

    test_thread_interning!(make_excitation_operator(0u32));

    // Implementation for Expr
    #[test]
    fn test_impl_expr() {
        let op1 = make_excitation_operator(0u32);
        let name = random_alphanumeric(0u32);

        let op = downcast_from_arc::<ExcitationOperator>(&op1).unwrap();
        assert_eq!(
            op,
            &ExcitationOperator {
                name: name.clone()
            }
        );
        assert_eq!(op.name(), name);

        assert_eq!(op1.hash_key(), format!("ExcitationOperator({name})"));
        assert!(!op1.is_scalar());
        assert_eq!(format!("{}", op1), name);

        let op2 = make_excitation_operator(0u32);
        let op3 = make_excitation_operator(3u32);
        assert_eq!(&op1, &op2);
        assert_ne!(&op1, &op3);
    }

    #[test]
    fn test_differentiation() {
        let op = make_excitation_operator(10u32);
        let p = make_perturbation_symbol(4u32, 4u32);
        assert_eq!(&op.differentiate(&p).unwrap(), &ZeroOperator::new());
    }

    // Test serialization and deserialization via `serde_json`
    #[test]
    fn test_serialization() {
        let op = make_excitation_operator(10u32);
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    // Test utils
    #[test]
    fn test_utils() {
        let op1 = make_excitation_operator(0u32);

        assert!(is_expr_type::<ExcitationOperator>(&op1));
        assert!(!is_zero_expr(&op1, None));
        assert!(!is_one_expr(&op1, None));

        let op2 = make_excitation_operator(0u32);
        let op3 = make_excitation_operator(10u32);

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));

        assert!(is_zero_expr(
            &MatrixAdd::new(vec![op1.clone(), negate_expr(op1.clone()).unwrap()]).unwrap(),
            None
        ));
        assert!(is_zero_expr(&subtract_exprs(op1.clone(), op2.clone()).unwrap(), None));
    }
}
