use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use typetag;

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::expressions::{LagMultiplier, WfnParameter, ZeroOperator};
use crate::internal::{intern_expr, multi_perturbation_format, multi_perturbation_hash};
use crate::perturbations::Perturbation;
use crate::public::{
    downcast_from_arc, downcast_from_ref, expression_error, generic_expression_error, is_expr_type,
};

/// A ResidueParameter is a perturbed parameter with the sum of frequencies of
/// some perturbations approaching the energy of an excited state from the
/// positive or negative side.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct ResidueParameter {
    positive_frequency: bool,
    perturbations: Vec<Arc<Perturbation>>,
    excited_state: Arc<dyn Expr>,
    parameter: Arc<dyn Expr>,
}

impl ResidueParameter {
    #[inline]
    pub fn builder(
        perturbations: Vec<Arc<Perturbation>>,
        excited_state: Arc<dyn Expr>,
        parameter: Arc<dyn Expr>,
    ) -> ResidueParameterBuilder {
        ResidueParameterBuilder {
            positive_frequency: true,
            perturbations,
            excited_state,
            parameter,
        }
    }

    #[inline]
    pub fn positive_frequency(&self) -> bool {
        self.positive_frequency
    }

    #[inline]
    pub fn perturbations(&self) -> &[Arc<Perturbation>] {
        &self.perturbations
    }

    #[inline]
    pub fn excited_state(&self) -> &Arc<dyn Expr> {
        &self.excited_state
    }

    #[inline]
    pub fn parameter(&self) -> &Arc<dyn Expr> {
        &self.parameter
    }
}

#[derive(Debug)]
pub struct ResidueParameterBuilder {
    positive_frequency: bool,
    perturbations: Vec<Arc<Perturbation>>,
    excited_state: Arc<dyn Expr>,
    parameter: Arc<dyn Expr>,
}

impl ResidueParameterBuilder {
    #[inline]
    pub fn positive_frequency(mut self, positive_frequency: bool) -> Self {
        self.positive_frequency = positive_frequency;
        self
    }

    pub fn build(mut self) -> Result<Arc<dyn Expr>, TinnedError> {
        if is_expr_type::<ZeroOperator>(&self.parameter) {
            return Ok(ZeroOperator::new());
        }

        let derivative_opt = if let Some(wfn) = downcast_from_arc::<WfnParameter>(&self.parameter) {
            Some(wfn.derivative())
        } else if let Some(lag) = downcast_from_arc::<LagMultiplier>(&self.parameter) {
            Some(lag.derivative())
        } else {
            None
        };

        let derivative = derivative_opt.ok_or_else(|| {
            expression_error(
                "ResidueParameterBuilder: parameter must be a WfnParameter or LagMultiplier",
                &self.parameter,
                None,
            )
        })?;

        self.perturbations.sort();

        if !derivative.is_subchain_vec(&self.perturbations) {
            return Ok(ZeroOperator::new());
        }

        Ok(intern_expr(Arc::new(ResidueParameter {
            positive_frequency: self.positive_frequency,
            perturbations: self.perturbations,
            excited_state: self.excited_state,
            parameter: self.parameter,
        })))
    }
}

impl ExprInternal for ResidueParameter {
    impl_expr_internal_methods!(ResidueParameter);

    #[inline]
    fn hash_key(&self) -> String {
        format!(
            "ResidueParameter([{}]; {}; {}; {})",
            multi_perturbation_hash(&self.perturbations, ";"),
            self.positive_frequency,
            self.excited_state.hash_key(),
            self.parameter.hash_key(),
        )
    }

    #[inline]
    fn total_order(&self) -> u32 {
        self.parameter.total_order()
    }

    #[inline]
    fn deep_eq_superchains(&self, other: &Arc<dyn Expr>) -> bool {
        if let Some(op) = downcast_from_arc::<ResidueParameter>(other) {
            self.parameter.deep_eq_superchains(&op.parameter)
        } else {
            false
        }
    }
}

#[typetag::serde]
impl Expr for ResidueParameter {
    // We treat `ResidueParameter` is the same type as its `parameter` so that
    // `exist_any()` will return true and `find_superchains()` will return
    // `ResidueParameter` itself if its `parameter` is the input parameter of
    // these two methods.
    //
    // Elimination should be usually performed for response functions, so that
    // it is a bit weird to call eliminate() method for `ResidueParameter`.
    //
    // But, it is technically possible to firstly take the residue procedure,
    // followed by the elimination in SymResponse. So we may still need the
    // eliminate() method for `ResidueParameter`.
    impl_unary_expr_common_methods!(
        ResidueParameter,
        parameter,
        False,
        false,
        |this: &ResidueParameter, arg| Self::builder(
            this.perturbations.clone(),
            this.excited_state.clone(),
            arg
        )
        .positive_frequency(this.positive_frequency)
        .build()
    );

    fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError> {
        let diff_param = self.parameter.differentiate(s).map_err(|e| {
            generic_expression_error(
                "ResidueParameter::differentiate() failed for parameter",
                self,
                Some(Box::new(e)),
            )
        })?;

        Self::builder(self.perturbations.clone(), self.excited_state.clone(), diff_param)
            .positive_frequency(self.positive_frequency)
            .build()
    }
}

impl PartialEq for ResidueParameter {
    fn eq(&self, other: &Self) -> bool {
        self.positive_frequency == other.positive_frequency
            && self.perturbations == other.perturbations
            && &self.excited_state == &other.excited_state
            && &self.parameter == &other.parameter
    }
}

impl Eq for ResidueParameter {}

impl std::fmt::Display for ResidueParameter {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "lim([{}]; {}{})({})",
            multi_perturbation_format(&self.perturbations, ";"),
            if self.positive_frequency {
                "-"
            } else {
                "+"
            },
            self.excited_state,
            self.parameter,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expressions::symbol::test_utils::make_symbol;
    //use crate::expressions::wfn_parameter::test_utils::make_wfn_parameter;
    use crate::perturbations::pert_multichain::test_utils::{make_pert_multichain, make_pert_vec};
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::{is_one_expr, is_zero_expr};

    test_struct_safety!(ResidueParameter);

    test_thread_interning!({
        let deriv = make_pert_multichain(0u32, 0u32, 2u32, 0u32);
        let perturbations = deriv.keys();
        let excited_state = make_symbol(0u32);
        let parameter = WfnParameter::builder("wfn").derivative(deriv).build().unwrap();
        ResidueParameter::builder(perturbations, excited_state, parameter)
            .positive_frequency(true)
            .build()
            .unwrap()
    });

    #[test]
    fn test_impl_expr() {
        let len_pert_name = 2u32;
        let deriv = make_pert_multichain(len_pert_name, 8u32, 1u32, 10u32);
        let mut perturbations = make_pert_vec(len_pert_name + 1, 8u32);
        let len_excited_name = 4u32;
        let excited_state = make_symbol(len_excited_name);
        let mut parameter = WfnParameter::builder("wfn").derivative(deriv.clone()).build().unwrap();
        let positive_frequency = true;
        let op1 = ResidueParameter::builder(
            perturbations.clone(),
            excited_state.clone(),
            parameter.clone(),
        )
        .build()
        .unwrap();

        assert!(is_expr_type::<ZeroOperator>(&op1));

        perturbations = deriv.keys();
        let op2 = ResidueParameter::builder(
            perturbations.clone(),
            excited_state.clone(),
            parameter.clone(),
        )
        .positive_frequency(positive_frequency)
        .build()
        .unwrap();

        let op = downcast_from_arc::<ResidueParameter>(&op2).unwrap();
        assert_eq!(
            op,
            &ResidueParameter {
                positive_frequency,
                perturbations: perturbations.clone(),
                excited_state: excited_state.clone(),
                parameter: parameter.clone(),
            }
        );

        assert_eq!(op.positive_frequency(), positive_frequency);
        assert_eq!(op.perturbations(), &perturbations);
        assert_eq!(op.excited_state(), &excited_state);
        assert_eq!(op.parameter(), &parameter);

        assert_eq!(
            op2.hash_key(),
            format!(
                "ResidueParameter([{}]; {}; {}; {})",
                multi_perturbation_hash(&perturbations, ";"),
                positive_frequency,
                excited_state.hash_key(),
                parameter.hash_key()
            )
        );
        assert!(!op2.is_scalar());
        assert_eq!(
            format!("{}", op2),
            format!(
                "lim([{}]; {}{})({})",
                multi_perturbation_format(&perturbations, ";"),
                if positive_frequency {
                    "-"
                } else {
                    "+"
                },
                excited_state,
                parameter,
            )
        );

        let op3 = ResidueParameter::builder(
            perturbations.clone(),
            excited_state.clone(),
            parameter.clone(),
        )
        .positive_frequency(positive_frequency)
        .build()
        .unwrap();

        assert_eq!(&op2, &op3);

        let op4 = ResidueParameter::builder(
            perturbations.clone(),
            excited_state.clone(),
            parameter.clone(),
        )
        .positive_frequency(!positive_frequency)
        .build()
        .unwrap();
        let op5 = ResidueParameter::builder(
            perturbations.clone(),
            make_symbol(len_excited_name + 1),
            parameter.clone(),
        )
        .positive_frequency(positive_frequency)
        .build()
        .unwrap();
        let p = make_perturbation_symbol(len_pert_name + 1, 8u32);
        parameter = parameter.differentiate(&p).unwrap();
        let op6 = ResidueParameter::builder(
            perturbations.clone(),
            excited_state.clone(),
            parameter.clone(),
        )
        .positive_frequency(positive_frequency)
        .build()
        .unwrap();
        perturbations.push(p.clone());
        let op7 = ResidueParameter::builder(
            perturbations.clone(),
            excited_state.clone(),
            parameter.clone(),
        )
        .positive_frequency(positive_frequency)
        .build()
        .unwrap();

        assert_ne!(&op2, &op4);
        assert_ne!(&op2, &op5);
        assert_ne!(&op2, &op6);
        assert_ne!(&op2, &op7);
    }

    #[test]
    fn test_differentiation() {
        let len_pert_name = 2u32;
        let deriv = make_pert_multichain(len_pert_name, 8u32, 1u32, 10u32);
        let parameter = WfnParameter::builder("wfn").derivative(deriv.clone()).build().unwrap();
        let op = ResidueParameter::builder(deriv.keys(), make_symbol(4u32), parameter.clone())
            .build()
            .unwrap();

        let p = make_perturbation_symbol(len_pert_name, 8u32);
        let diff_op = op.differentiate(&p).unwrap();
        let diff_param = parameter.differentiate(&p).unwrap();

        let op_cast = downcast_from_arc::<ResidueParameter>(&op).unwrap();
        let diff_cast = downcast_from_arc::<ResidueParameter>(&diff_op).unwrap();

        assert_eq!(diff_cast.positive_frequency(), op_cast.positive_frequency());
        assert_eq!(diff_cast.perturbations(), op_cast.perturbations());
        assert_eq!(diff_cast.excited_state(), op_cast.excited_state());
        assert_ne!(diff_cast.parameter(), op_cast.parameter());
        assert_eq!(diff_cast.parameter(), &diff_param);
    }

    #[test]
    fn test_serialization() {
        let deriv = make_pert_multichain(2u32, 8u32, 1u32, 10u32);
        let perturbations = deriv.keys();
        let excited_state = make_symbol(4u32);
        let parameter = WfnParameter::builder("wfn").derivative(deriv.clone()).build().unwrap();
        let positive_frequency = true;
        let op = ResidueParameter::builder(
            perturbations.clone(),
            excited_state.clone(),
            parameter.clone(),
        )
        .positive_frequency(positive_frequency)
        .build()
        .unwrap();
        let json = serde_json::to_string(&op).unwrap();
        let deserialized: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&op, &deserialized);
    }

    #[test]
    fn test_utils() {
        let len_pert_name = 2u32;
        let deriv = make_pert_multichain(len_pert_name, 8u32, 1u32, 10u32);
        let mut perturbations = deriv.keys();
        let len_excited_name = 4u32;
        let excited_state = make_symbol(len_excited_name);
        let mut parameter = WfnParameter::builder("wfn").derivative(deriv.clone()).build().unwrap();
        let positive_frequency = true;
        let op1 = ResidueParameter::builder(
            perturbations.clone(),
            excited_state.clone(),
            parameter.clone(),
        )
        .positive_frequency(positive_frequency)
        .build()
        .unwrap();

        assert!(is_expr_type::<ResidueParameter>(&op1));
        assert!(!is_zero_expr(&op1, None));
        assert!(!is_one_expr(&op1, None));

        let op2 = ResidueParameter::builder(
            perturbations.clone(),
            excited_state.clone(),
            parameter.clone(),
        )
        .positive_frequency(positive_frequency)
        .build()
        .unwrap();

        assert!(Arc::ptr_eq(&op1, &op2));

        let op4 = ResidueParameter::builder(
            perturbations.clone(),
            excited_state.clone(),
            parameter.clone(),
        )
        .positive_frequency(!positive_frequency)
        .build()
        .unwrap();
        let op5 = ResidueParameter::builder(
            perturbations.clone(),
            make_symbol(len_excited_name + 1),
            parameter.clone(),
        )
        .positive_frequency(positive_frequency)
        .build()
        .unwrap();
        let p = make_perturbation_symbol(len_pert_name + 1, 8u32);
        parameter = parameter.differentiate(&p).unwrap();
        let op6 = ResidueParameter::builder(
            perturbations.clone(),
            excited_state.clone(),
            parameter.clone(),
        )
        .positive_frequency(positive_frequency)
        .build()
        .unwrap();
        perturbations.push(p.clone());
        let op7 = ResidueParameter::builder(
            perturbations.clone(),
            excited_state.clone(),
            parameter.clone(),
        )
        .positive_frequency(positive_frequency)
        .build()
        .unwrap();

        assert!(!Arc::ptr_eq(&op2, &op4));
        assert!(!Arc::ptr_eq(&op2, &op5));
        assert!(!Arc::ptr_eq(&op2, &op6));
        assert!(!Arc::ptr_eq(&op2, &op7));
    }
}
