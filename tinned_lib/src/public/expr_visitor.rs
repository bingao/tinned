// This is a direct benefit of Rust's monomorphization and lazy compilation --
// unused generic code does not cause a build failure.

// Clearly document that accept() is a generic hook -- users can ignore it
// unless they need evaluation.
//
// Consider adding a test evaluator in your own crate to ensure correctness of
// visitor plumbing, but behind a dev-only feature flag, so it doesn't bloat
// your public API.

pub trait ExprVisitor {
    type Output;

    fn visit_number(&mut self, n: &Number) -> Result<Self::Output, EvalError>;
    fn visit_add(&mut self, add: &Add) -> Result<Self::Output, EvalError>;
    fn visit_mul(&mut self, mul: &Mul) -> Result<Self::Output, EvalError>;
    // Add other visit_* methods as needed
}

// Step 2: Add accept() to the Expr trait
pub trait Expr: Debug + Send + Sync {
    fn accept<V: crate::public::expr_visitor::ExprVisitor>(&self, visitor: &mut V) -> Result<V::Output, EvalError>;
}

// Step 3: Implement accept in each concrete expression type
impl Expr for Add {
    fn accept<V: ExprVisitor>(&self, visitor: &mut V) -> Result<V::Output, EvalError> {
        visitor.visit_add(self)
    }
}

impl Expr for Number {
    fn accept<V: ExprVisitor>(&self, visitor: &mut V) -> Result<V::Output, EvalError> {
        visitor.visit_number(self)
    }
}

// Re-export for user convenience
// In src/lib.rs: pub use public::expr_visitor::ExprVisitor;

// Step 4: Users implement their own evaluator

pub struct NumericEvaluator;

impl ExprVisitor for NumericEvaluator {
    type Output = f64;

    fn visit_number(&mut self, n: &Number) -> Result<Self::Output, EvalError> {
        n.as_f64() // for example
    }

    fn visit_add(&mut self, add: &Add) -> Result<Self::Output, EvalError> {
        add.terms()
            .iter()
            .map(|term| term.accept(self))
            .sum()
    }

    fn visit_mul(&mut self, mul: &Mul) -> Result<Self::Output, EvalError> {
        mul.terms()
            .iter()
            .map(|term| term.accept(self))
            .product()
    }
}

// Step 5: Evaluation entry point

let expr: Arc<dyn Expr> = ...;

let mut evaluator = NumericEvaluator;
let result: f64 = expr.accept(&mut evaluator)?;

