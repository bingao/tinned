use std::sync::Arc;

use crate::core::{Expr, TinnedError};
use crate::expressions::{
    Add, AdjointMap, Composition, Conjugate, DotProduct, ExchCorrEnergy, ExchCorrPotential,
    ExpAdjointMap, HermitianTranspose, LagMultiplier, MatrixAdd, MatrixMul, Mul, NonElecFunction,
    Number, OneElecOperator, Power, ResidueParameter, Symbol, TemporumOperator, TemporumOverlap,
    Trace, Transpose, TwoElecEnergy, TwoElecOperator, WfnParameter, ZeroOperator,
};
use crate::public::downcast_from_arc;

#[cfg_attr(feature = "ffi", derive(safer_ffi::derive_ReprC))]
#[cfg_attr(feature = "ffi", repr(C))]
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ExprTag {
    Add,
    AdjointMap,
    Composition,
    Conjugate,
    DotProduct,
    ExchCorrEnergy,
    ExchCorrPotential,
    ExpAdjointMap,
    HermitianTranspose,
    LagMultiplier,
    MatrixAdd,
    MatrixMul,
    Mul,
    NonElecFunction,
    Number,
    OneElecOperator,
    Power,
    ResidueParameter,
    Symbol,
    TemporumOperator,
    TemporumOverlap,
    Trace,
    Transpose,
    TwoElecEnergy,
    TwoElecOperator,
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
        if downcast_from_arc::<Number>(expr).is_some() {
            return visitor.leaf(Number, expr);
        }
        if downcast_from_arc::<Symbol>(expr).is_some() {
            return visitor.leaf(Symbol, expr);
        }
        if downcast_from_arc::<ZeroOperator>(expr).is_some() {
            return visitor.leaf(ZeroOperator, expr);
        }

        // Composite expressions
        if let Some(add) = downcast_from_arc::<Add>(expr) {
            visitor.begin(Add, add.terms().len())?;
            for term in add.terms() {
                walk(term, visitor)?;
            }
            return visitor.end(Add, add.terms().len());
        }
        if let Some(mul) = downcast_from_arc::<Mul>(expr) {
            visitor.begin(Mul, mul.factors().len()+1)?;
            walk(mul.coefficient().into(), visitor)?;
            for factor in mul.factors() {
                walk(factor, visitor)?;
            }
            return visitor.end(Mul, mul.factors().len()+1);
        }
        //if let Some(p) = downcast_from_arc::<Power>(e) {
        //    v.begin(Power, 2)?;
        //    walk(p.base(), v)?; // base
        //    walk(p.exponent(), v)?; // exponent (adjust if exponent is int)
        //    return v.end(Power, 2);
        //}

        Ok(())
    }

    walk(root, visitor)
}
