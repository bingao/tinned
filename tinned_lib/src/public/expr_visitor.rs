use std::sync::Arc;

use crate::core::{Expr, TinnedError};
use crate::expressions::{
    Add, AdjointMap, AoTwoElecEnergy, AoTwoElecMatrix, BasisTimeEvolution, Composition, Conjugate,
    DotProduct, ExchCorrEnergy, ExchCorrPotential, ExcitationOperator, ExpAdjointMap,
    LagMultiplier, MatrixAdd, MatrixMul, Mul, NonElecFunction, Number, OneElecMatrix, Power,
    ResidueParameter, SubExpr, Symbol, TimeEvolution, Trace, Transpose, TwoElecMatrix,
    WfnParameter, ZeroOperator,
};
use crate::public::downcast_from_arc;

#[cfg_attr(feature = "ffi", safer_ffi::derive_ReprC)]
#[repr(u32)]
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ExprTag {
    Add,
    AdjointMap,
    AoTwoElecEnergy,
    AoTwoElecMatrix,
    BasisTimeEvolution,
    Composition,
    Conjugate,
    DotProduct,
    ExchCorrEnergy,
    ExchCorrPotential,
    ExcitationOperator,
    ExpAdjointMap,
    LagMultiplier,
    MatrixAdd,
    MatrixMul,
    Mul,
    NonElecFunction,
    Number,
    OneElecMatrix,
    Power,
    ResidueParameter,
    SubExpr,
    Symbol,
    TimeEvolution,
    Trace,
    Transpose,
    TwoElecMatrix,
    WfnParameter,
    ZeroOperator,
}

// A general expression visitor that users need to develop their own ones
pub trait ExprVisitor {
    // Called for expressions that have children to iterate, such as addition and multiplication
    fn begin(&mut self, tag: ExprTag, arity: usize) -> Result<(), TinnedError>;

    // Called for leaf expressions, which users can access their fields and perform corresponding actions
    fn leaf(&mut self, tag: ExprTag, expr: &Arc<dyn Expr>) -> Result<(), TinnedError>;

    // Called after all children being visited
    fn end(&mut self, tag: ExprTag, arity: usize) -> Result<(), TinnedError>;
}

// Visit an expression using post-order traversal
pub fn walk_expr_postorder<V: ExprVisitor>(
    root: &Arc<dyn Expr>,
    visitor: &mut V,
) -> Result<(), TinnedError> {
    fn walk<V: ExprVisitor>(expr: &Arc<dyn Expr>, visitor: &mut V) -> Result<(), TinnedError> {
        // Leaf expressions
        if downcast_from_arc::<AdjointMap>(expr).is_some() {
            return visitor.leaf(ExprTag::AdjointMap, expr);
        }
        if downcast_from_arc::<AoTwoElecEnergy>(expr).is_some() {
            return visitor.leaf(ExprTag::AoTwoElecEnergy, expr);
        }
        if downcast_from_arc::<AoTwoElecMatrix>(expr).is_some() {
            return visitor.leaf(ExprTag::AoTwoElecMatrix, expr);
        }
        if downcast_from_arc::<BasisTimeEvolution>(expr).is_some() {
            return visitor.leaf(ExprTag::BasisTimeEvolution, expr);
        }
        if downcast_from_arc::<Composition>(expr).is_some() {
            return visitor.leaf(ExprTag::Composition, expr);
        }
        if downcast_from_arc::<Conjugate>(expr).is_some() {
            return visitor.leaf(ExprTag::Conjugate, expr);
        }
        if downcast_from_arc::<DotProduct>(expr).is_some() {
            return visitor.leaf(ExprTag::DotProduct, expr);
        }
        if downcast_from_arc::<ExchCorrEnergy>(expr).is_some() {
            return visitor.leaf(ExprTag::ExchCorrEnergy, expr);
        }
        if downcast_from_arc::<ExchCorrPotential>(expr).is_some() {
            return visitor.leaf(ExprTag::ExchCorrPotential, expr);
        }
        if downcast_from_arc::<ExcitationOperator>(expr).is_some() {
            return visitor.leaf(ExprTag::ExcitationOperator, expr);
        }
        if downcast_from_arc::<ExpAdjointMap>(expr).is_some() {
            return visitor.leaf(ExprTag::ExpAdjointMap, expr);
        }
        if downcast_from_arc::<LagMultiplier>(expr).is_some() {
            return visitor.leaf(ExprTag::LagMultiplier, expr);
        }
        if downcast_from_arc::<NonElecFunction>(expr).is_some() {
            return visitor.leaf(ExprTag::NonElecFunction, expr);
        }
        if downcast_from_arc::<Number>(expr).is_some() {
            return visitor.leaf(ExprTag::Number, expr);
        }
        if downcast_from_arc::<OneElecMatrix>(expr).is_some() {
            return visitor.leaf(ExprTag::OneElecMatrix, expr);
        }
        if downcast_from_arc::<Power>(expr).is_some() {
            return visitor.leaf(ExprTag::Power, expr);
        }
        if downcast_from_arc::<ResidueParameter>(expr).is_some() {
            return visitor.leaf(ExprTag::ResidueParameter, expr);
        }
        if downcast_from_arc::<Symbol>(expr).is_some() {
            return visitor.leaf(ExprTag::Symbol, expr);
        }
        if downcast_from_arc::<TimeEvolution>(expr).is_some() {
            return visitor.leaf(ExprTag::TimeEvolution, expr);
        }
        if downcast_from_arc::<Trace>(expr).is_some() {
            return visitor.leaf(ExprTag::Trace, expr);
        }
        if downcast_from_arc::<Transpose>(expr).is_some() {
            return visitor.leaf(ExprTag::Transpose, expr);
        }
        if downcast_from_arc::<TwoElecMatrix>(expr).is_some() {
            return visitor.leaf(ExprTag::TwoElecMatrix, expr);
        }
        if downcast_from_arc::<WfnParameter>(expr).is_some() {
            return visitor.leaf(ExprTag::WfnParameter, expr);
        }
        if downcast_from_arc::<ZeroOperator>(expr).is_some() {
            return visitor.leaf(ExprTag::ZeroOperator, expr);
        }

        // Composite expressions
        if let Some(add) = downcast_from_arc::<Add>(expr) {
            visitor.begin(ExprTag::Add, add.terms().len())?;
            for term in add.terms() {
                walk(term, visitor)?;
            }
            return visitor.end(ExprTag::Add, add.terms().len());
        }
        if let Some(mul) = downcast_from_arc::<Mul>(expr) {
            visitor.begin(ExprTag::Mul, mul.factors().len() + 1)?;
            let coefficient: Arc<dyn Expr> = mul.coefficient().into();
            walk(&coefficient, visitor)?;
            for factor in mul.factors() {
                walk(factor, visitor)?;
            }
            return visitor.end(ExprTag::Mul, mul.factors().len() + 1);
        }
        if let Some(add) = downcast_from_arc::<MatrixAdd>(expr) {
            visitor.begin(ExprTag::MatrixAdd, add.terms().len())?;
            for term in add.terms() {
                walk(term, visitor)?;
            }
            return visitor.end(ExprTag::MatrixAdd, add.terms().len());
        }
        if let Some(mul) = downcast_from_arc::<MatrixMul>(expr) {
            visitor.begin(ExprTag::MatrixMul, mul.factors().len() + 1)?;
            walk(mul.coefficient(), visitor)?;
            for factor in mul.factors() {
                walk(factor, visitor)?;
            }
            return visitor.end(ExprTag::MatrixMul, mul.factors().len() + 1);
        }
        if let Some(sub_expr) = downcast_from_arc::<SubExpr>(expr) {
            visitor.begin(ExprTag::SubExpr, 1)?;
            walk(sub_expr.expression(), visitor)?;
            return visitor.end(ExprTag::Add, 1);
        }

        Ok(())
    }

    walk(root, visitor)
}
