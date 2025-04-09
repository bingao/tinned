use std::any::Any;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::ops::{Add as StdAdd, Mul as StdMul};
use std::sync::Arc;

use num_complex::Complex64;
use num_rational::Rational64;
use serde::{Deserialize, Serialize};

use crate::core::{Expr, TinnedError};
use crate::perturbations::Perturbation;
use crate::utils::intern;

// Define an enum to store different number types
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Number {
    Integer(i64),
    Real(f64),
    Complex(Complex64),
    Fraction(Rational64),
}

impl Number {
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

// From implementations for primitives, so one can use, for example,
//
// let a: Arc<dyn Expr> = 3.into();
// let b: Arc<dyn Expr> = 2.0.into();
// let c: Arc<dyn Expr> = Complex64::new(1.0, 2.0).into();
macro_rules! impl_from_number {
    ($t:ty, $variant:ident) => {
        impl From<$t> for Arc<dyn Expr> {
            #[inline]
            fn from(value: $t) -> Self {
                intern(Arc::new(Number::$variant(value)))
            }
        }
    };
}

impl_from_number!(i64, Integer);
impl_from_number!(f64, Real);
impl_from_number!(Complex64, Complex);
impl_from_number!(Rational64, Fraction);

// From<Number> implementation
impl From<Number> for Arc<dyn Expr> {
    #[inline]
    fn from(num: Number) -> Self {
        intern(Arc::new(num))
    }
}

impl Expr for Number {
    #[inline]
    fn as_any(&self) -> &dyn Any {
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
    fn differentiate(&self, _s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(0.into())
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

impl Display for Number {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        match self {
            Number::Integer(n) => write!(f, "{}", n),
            Number::Real(r) => write!(f, "{}", r),
            Number::Complex(c) => write!(f, "{} + {}i", c.re, c.im),
            Number::Fraction(r) => write!(f, "{}/{}", r.numer(), r.denom()),
        }
    }
}
