use std::sync::Arc;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::expressions::{OneElecOperator, WfnParameter, ZeroOperator};
use crate::utils::{
    downcast_from_ref, intern_expr, invalid_expression_error, is_expr_type, is_zero_expr,
};

/// A TemporumOperator is a non-scalar operator acting on a ket or a bra
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct TemporumOperator {
    on_ket: bool,
    argument: Arc<dyn Expr>,
}

impl TemporumOperator {
    #[inline]
    pub fn builder(argument: Arc<dyn Expr>) -> TemporumOperatorBuilder {
        TemporumOperatorBuilder {
            on_ket: true,
            argument,
        }
    }

    #[inline]
    fn builder_from(&self, argument: Arc<dyn Expr>) -> TemporumOperatorBuilder {
        TemporumOperatorBuilder {
            on_ket: self.on_ket,
            argument,
        }
    }

    #[inline]
    pub fn on_ket(&self) -> bool {
        self.on_ket
    }

    #[inline]
    pub fn argument(&self) -> &Arc<dyn Expr> {
        &self.argument
    }
}

#[derive(Debug)]
pub struct TemporumOperatorBuilder {
    on_ket: bool,
    argument: Arc<dyn Expr>,
}

impl TemporumOperatorBuilder {
    #[inline]
    pub fn on_ket(mut self, on_ket: bool) -> Self {
        self.on_ket = on_ket;
        self
    }

    pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.argument.is_scalar() {
            return Err(invalid_expression_error(
                "TemporumOperatorBuilder::build() - scalar argument",
                &self.argument,
            ));
        }

        if is_expr_type::<OneElecOperator>(&self.argument)
            || is_expr_type::<WfnParameter>(&self.argument)
        {
            Ok(intern_expr(Arc::new(TemporumOperator {
                on_ket: self.on_ket,
                argument: self.argument,
            })))
        } else {
            Err(invalid_expression_error(
                "TemporumOperatorBuilder::build() - unsupported argument type",
                &self.argument,
            ))
        }
    }
}

#[typetag::serde]
impl Expr for TemporumOperator {
    #[inline]
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        format!("TemporumOperator({}; {})", self.on_ket, self.argument.hash_key())
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        false
    }

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        if let Some(op) = downcast_from_ref::<TemporumOperator>(other) {
            self == op
        } else {
            false
        }
    }

    #[inline]
    fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }

    fn differentiate(
        &self,
        s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_arg = self.argument.differentiate(s)?;

        if is_zero_expr(&diff_arg) {
            Ok(ZeroOperator::new())
        } else {
            self.builder_from(diff_arg).build()
        }
    }
}

impl PartialEq for TemporumOperator {
    fn eq(&self, other: &Self) -> bool {
        self.on_ket == other.on_ket && &self.argument == &other.argument
    }
}

impl Eq for TemporumOperator {}

impl std::fmt::Display for TemporumOperator {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}({})",
            if self.on_ket {
                "i*dt"
            } else {
                "-i*dt"
            },
            self.argument
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::one_elec_operator::test_utils::make_one_elec_operator;
    use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::utils::{downcast_from_arc, is_one_expr};

    test_struct_safety!(TemporumOperator);

    test_thread_interning!(TemporumOperator::builder(make_one_elec_operator("1el"))
        .build()
        .unwrap());

    #[test]
    fn test_impl_expr() {
        let on_ket = true;
        let argument = make_one_elec_operator("");
        let op1 = TemporumOperator::builder(argument.clone()).on_ket(on_ket).build().unwrap();

        let op = downcast_from_arc::<TemporumOperator>(&op1).unwrap();
        assert_eq!(
            op,
            &TemporumOperator {
                on_ket,
                argument: argument.clone()
            }
        );

        assert_eq!(op.on_ket(), on_ket);
        assert_eq!(op.argument(), &argument.clone());

        let op2 = op.builder_from(argument.clone()).build().unwrap();
        assert!(Arc::ptr_eq(&op1, &op2));
        assert_eq!(&op1, &op2);

        assert_eq!(
            op1.hash_key(),
            format!("TemporumOperator({}; {})", on_ket, argument.hash_key())
        );
        assert!(!op1.is_scalar());
        assert_eq!(
            format!("{}", op1),
            format!(
                "{}({})",
                if on_ket {
                    "i*dt"
                } else {
                    "-i*dt"
                },
                argument
            )
        );

        let op3 = TemporumOperator::builder(argument).on_ket(!on_ket).build().unwrap();
        let op4 = TemporumOperator::builder(make_wfn_parameter("")).build().unwrap();

        assert_ne!(&op1, &op3);
        assert_ne!(&op1, &op4);
    }

    #[test]
    fn test_serialization() {
        let mut op = TemporumOperator::builder(make_one_elec_operator("")).build().unwrap();
        let mut json = serde_json::to_string(&op).unwrap();
        let mut deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);

        op = TemporumOperator::builder(make_wfn_parameter("")).build().unwrap();
        json = serde_json::to_string(&op).unwrap();
        deserialized = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let op1 = TemporumOperator::builder(make_one_elec_operator("1el")).build().unwrap();

        assert!(is_expr_type::<TemporumOperator>(&op1));
        assert!(!is_zero_expr(&op1));
        assert!(!is_one_expr(&op1));

        let op2 = TemporumOperator::builder(make_one_elec_operator("1el")).build().unwrap();
        let op3 = TemporumOperator::builder(make_one_elec_operator("")).build().unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));
        assert!(!Arc::ptr_eq(&op1, &op3));
    }
}
