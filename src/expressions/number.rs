use std::sync::Arc;

use num_complex::Complex64;
use num_rational::Rational64;

use typetag;

use crate::core::{Expr, TinnedError};
use crate::utils::intern_expr;

// Define an enum to store different number types
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub enum Number {
    Integer(i64),
    Real(f64),
    Complex(Complex64),
    Fraction(Rational64),
}

impl Number {
    #[inline]
    pub fn from_i64(value: i64) -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Integer(value)))
    }

    #[inline]
    pub fn from_f64(value: f64) -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Real(value)))
    }

    #[inline]
    pub fn from_complex(value: Complex64) -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Complex(value)))
    }

    #[inline]
    pub fn from_rational(value: Rational64) -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Fraction(value)))
    }

    #[inline]
    pub fn zero() -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Integer(0)))
    }

    #[inline]
    pub fn one() -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Integer(1)))
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        match self {
            Number::Integer(n) => *n == 0,
            Number::Real(n) => *n == 0.0,
            Number::Fraction(n) => *n == Rational64::from_integer(0),
            Number::Complex(c) => c.re == 0.0 && c.im == 0.0,
        }
    }

    #[inline]
    pub fn is_one(&self) -> bool {
        match self {
            Number::Integer(n) => *n == 1,
            Number::Real(n) => *n == 1.0,
            Number::Fraction(n) => *n == Rational64::from_integer(1),
            Number::Complex(c) => c.re == 1.0 && c.im == 0.0,
        }
    }

    #[inline]
    pub fn conjugate(&self) -> Number {
        match self {
            Number::Complex(z) => Number::Complex(z.conj()),
            _ => self.clone(),
        }
    }

    #[inline]
    pub fn add(&self, other: &Number) -> Number {
        use Number::*;

        match (self, other) {
            (Integer(a), Integer(b)) => Integer(a + b),

            (Integer(a), Fraction(b)) => Fraction(*b + Rational64::from(*a)),
            (Fraction(a), Integer(b)) => Fraction(*a + Rational64::from(*b)),
            (Fraction(a), Fraction(b)) => Fraction(*a + *b),

            (Integer(a), Real(b)) => Real(*a as f64 + *b),
            (Real(a), Integer(b)) => Real(*a + *b as f64),
            (Real(a), Real(b)) => Real(*a + *b),

            (Integer(a), Complex(b)) => Complex(Complex64::new(*a as f64, 0.0) + *b),
            (Complex(a), Integer(b)) => Complex(*a + Complex64::new(*b as f64, 0.0)),

            (Real(a), Fraction(b)) => Real(*a + b.to_f64()),
            (Fraction(a), Real(b)) => Real(a.to_f64() + *b),

            (Real(a), Complex(b)) => Complex(Complex64::new(*a, 0.0) + *b),
            (Complex(a), Real(b)) => Complex(*a + Complex64::new(*b, 0.0)),

            (Fraction(a), Complex(b)) => Complex(Complex64::new(a.to_f64(), 0.0) + *b),
            (Complex(a), Fraction(b)) => Complex(*a + Complex64::new(b.to_f64(), 0.0)),

            (Complex(a), Complex(b)) => Complex(*a + *b),
        }
    }

    #[inline]
    pub fn mul(&self, other: &Number) -> Number {
        use Number::*;

        match (self, other) {
            (Integer(a), Integer(b)) => Integer(a * b),

            (Integer(a), Fraction(b)) => Fraction(*b * Rational64::from(*a)),
            (Fraction(a), Integer(b)) => Fraction(*a * Rational64::from(*b)),
            (Fraction(a), Fraction(b)) => Fraction(*a * *b),

            (Integer(a), Real(b)) => Real(*a as f64 * *b),
            (Real(a), Integer(b)) => Real(*a * *b as f64),
            (Real(a), Real(b)) => Real(*a * *b),

            (Integer(a), Complex(b)) => Complex(Complex64::new(*a as f64, 0.0) * *b),
            (Complex(a), Integer(b)) => Complex(*a * Complex64::new(*b as f64, 0.0)),

            (Real(a), Fraction(b)) => Real(*a * b.to_f64()),
            (Fraction(a), Real(b)) => Real(a.to_f64() * *b),

            (Real(a), Complex(b)) => Complex(Complex64::new(*a, 0.0) * *b),
            (Complex(a), Real(b)) => Complex(*a * Complex64::new(*b, 0.0)),

            (Fraction(a), Complex(b)) => Complex(Complex64::new(a.to_f64(), 0.0) * *b),
            (Complex(a), Fraction(b)) => Complex(*a * Complex64::new(b.to_f64(), 0.0)),

            (Complex(a), Complex(b)) => Complex(*a * *b),
        }
    }
}

// From<Number> implementation so that we can use into() method for Number
impl From<Number> for Arc<dyn Expr> {
    #[inline]
    fn from(num: Number) -> Self {
        crate::utils::intern_expr(Arc::new(num))
    }
}

#[typetag::serde]
impl Expr for Number {
    #[inline]
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    #[inline]
    fn hash_key(&self) -> String {
        match self {
            Number::Integer(n) => format!("Integer({})", n),
            Number::Real(r) => format!("Real({})", r),
            Number::Complex(c) => format!("Complex({})", c),
            Number::Fraction(f) => format!("Fraction({}/{})", f.numer(), f.denom()),
        }
    }

    #[inline]
    fn is_scalar(&self) -> bool {
        true
    }

    #[inline]
    fn eq_expr(&self, other: &dyn Expr) -> bool {
        if let Some(num) = crate::utils::downcast_from_ref::<Number>(other) {
            self == num
        } else {
            false
        }
    }

    #[inline]
    fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }

    #[inline]
    fn differentiate(
        &self,
        _s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(Number::zero())
    }
}

impl PartialEq for Number {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Number::Integer(a), Number::Integer(b)) => a == b,
            (Number::Real(a), Number::Real(b)) => a == b,
            (Number::Complex(a), Number::Complex(b)) => a == b,
            (Number::Fraction(a), Number::Fraction(b)) => a == b,
            _ => false,
        }
    }
}

impl Eq for Number {}

impl std::fmt::Display for Number {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Number::Integer(n) => write!(f, "{}", n),
            Number::Real(r) => write!(f, "{}", r),
            Number::Complex(c) => write!(f, "{} + {}i", c.re, c.im),
            Number::Fraction(r) => write!(f, "{}/{}", r.numer(), r.denom()),
        }
    }
}
