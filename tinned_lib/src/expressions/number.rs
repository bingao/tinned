use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use float_cmp::approx_eq;
use num_complex::Complex64;
use num_rational::Rational64;
use num_traits::{ToPrimitive, Zero};

use crate::core::expr_internal::sealed::ExprInternal;
use crate::core::{Expr, TinnedError};
use crate::internal::intern_expr;
use crate::public::{NumberTolerance, generic_error, get_number_tolerance};

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
    pub fn from_i64(n: i64) -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Integer(n)))
    }

    #[inline]
    pub fn from_f64(f: f64) -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Real(f)))
    }

    #[inline]
    pub fn from_complex(z: Complex64) -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Complex(z)))
    }

    #[inline]
    pub fn from_rational(r: Rational64) -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Fraction(r)))
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
    pub fn minus_one() -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Integer(-1)))
    }

    #[inline]
    pub fn one_half() -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Fraction(Rational64::new(1, 2))))
    }

    #[inline]
    pub fn imaginary_unit() -> Arc<dyn Expr> {
        intern_expr(Arc::new(Number::Complex(Complex64::new(0.0, 1.0))))
    }

    #[inline]
    pub fn is_zero(&self, num_tol: Option<NumberTolerance>) -> bool {
        self.approx_eq_number(&Number::Integer(0), num_tol)
    }

    #[inline]
    pub fn is_one(&self, num_tol: Option<NumberTolerance>) -> bool {
        self.approx_eq_number(&Number::Integer(1), num_tol)
    }

    #[inline]
    pub fn conjugate(&self) -> Number {
        match self {
            Number::Complex(z) => Number::Complex(z.conj()),
            _ => self.clone(),
        }
    }

    #[inline]
    pub fn negate(&self) -> Number {
        match self {
            Number::Integer(n) => Number::Integer(-n),
            Number::Real(f) => Number::Real(-f),
            Number::Complex(z) => Number::Complex(Complex64::new(-z.re, -z.im)),
            Number::Fraction(r) => Number::Fraction(Rational64::new(-r.numer(), *r.denom())),
        }
    }

    #[inline]
    pub fn approx_eq_number(&self, other: &Number, num_tol: Option<NumberTolerance>) -> bool {
        use Number::*;

        let tol = num_tol.unwrap_or_else(|| get_number_tolerance());

        match (self, other) {
            (Integer(a), Integer(b)) => a == b,
            (Integer(a), Real(b)) => {
                approx_eq!(f64, *a as f64, *b, epsilon = tol.max_abs_error(*a as f64, *b))
            },
            (Integer(a), Complex(b)) => {
                approx_eq!(f64, *a as f64, b.re, epsilon = tol.max_abs_error(*a as f64, b.re))
                    && approx_eq!(f64, 0.0, b.im, epsilon = tol.max_abs_error(0.0, b.im))
            },
            (Integer(a), Fraction(b)) => Rational64::from_integer(*a) == *b,

            (Real(a), Integer(b)) => {
                approx_eq!(f64, *a, *b as f64, epsilon = tol.max_abs_error(*a, *b as f64))
            },
            (Real(a), Real(b)) => {
                approx_eq!(f64, *a, *b, epsilon = tol.max_abs_error(*a, *b))
            },
            (Real(a), Complex(b)) => {
                approx_eq!(f64, *a, b.re, epsilon = tol.max_abs_error(*a, b.re))
                    && approx_eq!(f64, 0.0, b.im, epsilon = tol.max_abs_error(0.0, b.im))
            },
            (Real(a), Fraction(b)) => {
                let b_f64 = b.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", b));
                approx_eq!(f64, *a, b_f64, epsilon = tol.max_abs_error(*a, b_f64))
            },

            (Complex(a), Integer(b)) => {
                approx_eq!(f64, a.re, *b as f64, epsilon = tol.max_abs_error(a.re, *b as f64))
                    && approx_eq!(f64, a.im, 0.0, epsilon = tol.max_abs_error(a.im, 0.0))
            },
            (Complex(a), Real(b)) => {
                approx_eq!(f64, a.re, *b, epsilon = tol.max_abs_error(a.re, *b))
                    && approx_eq!(f64, a.im, 0.0, epsilon = tol.max_abs_error(a.im, 0.0))
            },
            (Complex(a), Complex(b)) => {
                approx_eq!(f64, a.re, b.re, epsilon = tol.max_abs_error(a.re, b.re))
                    && approx_eq!(f64, a.im, b.im, epsilon = tol.max_abs_error(a.im, b.im))
            },
            (Complex(a), Fraction(b)) => {
                let b_f64 = b.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", b));
                approx_eq!(f64, a.re, b_f64, epsilon = tol.max_abs_error(a.re, b_f64))
                    && approx_eq!(f64, a.im, 0.0, epsilon = tol.max_abs_error(a.im, 0.0))
            },

            (Fraction(a), Integer(b)) => *a == Rational64::from_integer(*b),
            (Fraction(a), Real(b)) => {
                let a_f64 = a.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", a));
                approx_eq!(f64, a_f64, *b, epsilon = tol.max_abs_error(a_f64, *b))
            },
            (Fraction(a), Complex(b)) => {
                let a_f64 = a.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", a));
                approx_eq!(f64, a_f64, b.re, epsilon = tol.max_abs_error(a_f64, b.re))
                    && approx_eq!(f64, 0.0, b.im, epsilon = tol.max_abs_error(0.0, b.im))
            },
            (Fraction(a), Fraction(b)) => *a == *b,
        }
    }

    #[inline]
    pub fn add(&self, other: &Number) -> Number {
        use Number::*;

        match (self, other) {
            (Integer(a), Integer(b)) => Integer(a + b),
            (Integer(a), Real(b)) => Real(*a as f64 + *b),
            (Integer(a), Complex(b)) => Complex(Complex64::new(*a as f64, 0.0) + *b),
            (Integer(a), Fraction(b)) => Fraction(*b + Rational64::from_integer(*a)),

            (Real(a), Integer(b)) => Real(*a + *b as f64),
            (Real(a), Real(b)) => Real(*a + *b),
            (Real(a), Complex(b)) => Complex(Complex64::new(*a, 0.0) + *b),
            (Real(a), Fraction(b)) => {
                Real(*a + b.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", b)))
            },

            (Complex(a), Integer(b)) => Complex(*a + Complex64::new(*b as f64, 0.0)),
            (Complex(a), Real(b)) => Complex(*a + Complex64::new(*b, 0.0)),
            (Complex(a), Complex(b)) => Complex(*a + *b),
            (Complex(a), Fraction(b)) => Complex(
                *a + Complex64::new(
                    b.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", b)),
                    0.0,
                ),
            ),

            (Fraction(a), Integer(b)) => Fraction(*a + Rational64::from_integer(*b)),
            (Fraction(a), Real(b)) => {
                Real(a.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", a)) + *b)
            },
            (Fraction(a), Complex(b)) => Complex(
                Complex64::new(
                    a.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", a)),
                    0.0,
                ) + *b,
            ),
            (Fraction(a), Fraction(b)) => Fraction(*a + *b),
        }
    }

    #[inline]
    pub fn mul(&self, other: &Number) -> Number {
        use Number::*;

        match (self, other) {
            (Integer(a), Integer(b)) => Integer(a * b),
            (Integer(a), Real(b)) => Real(*a as f64 * *b),
            (Integer(a), Complex(b)) => Complex(Complex64::new(*a as f64, 0.0) * *b),
            (Integer(a), Fraction(b)) => Fraction(*b * Rational64::from_integer(*a)),

            (Real(a), Integer(b)) => Real(*a * *b as f64),
            (Real(a), Real(b)) => Real(*a * *b),
            (Real(a), Complex(b)) => Complex(Complex64::new(*a, 0.0) * *b),
            (Real(a), Fraction(b)) => {
                Real(*a * b.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", b)))
            },

            (Complex(a), Integer(b)) => Complex(*a * Complex64::new(*b as f64, 0.0)),
            (Complex(a), Real(b)) => Complex(*a * Complex64::new(*b, 0.0)),
            (Complex(a), Complex(b)) => Complex(*a * *b),
            (Complex(a), Fraction(b)) => Complex(
                *a * Complex64::new(
                    b.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", b)),
                    0.0,
                ),
            ),

            (Fraction(a), Integer(b)) => Fraction(*a * Rational64::from_integer(*b)),
            (Fraction(a), Real(b)) => {
                Real(a.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", a)) * *b)
            },
            (Fraction(a), Complex(b)) => Complex(
                Complex64::new(
                    a.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", a)),
                    0.0,
                ) * *b,
            ),
            (Fraction(a), Fraction(b)) => Fraction(*a * *b),
        }
    }

    #[inline]
    pub fn pow_i64(&self, exp: i64) -> Result<Number, TinnedError> {
        use Number::*;

        match self {
            Integer(n) => {
                if exp >= 0 {
                    Ok(Integer(n.pow(exp as u32)))
                } else {
                    // Negative power: promote to Fraction
                    if *n == 0 {
                        Err(generic_error("Cannot raise zero integer to negative power", None))
                    } else {
                        Ok(Fraction(
                            Rational64::from_integer(1)
                                / Rational64::from_integer(n.pow((-exp) as u32)),
                        ))
                    }
                }
            },

            Real(f) => {
                if *f == 0.0 && exp < 0 {
                    Err(generic_error("Cannot raise zero real number to negative power", None))
                } else {
                    Ok(Real(f.powi(exp as i32)))
                }
            },

            Complex(z) => {
                if z.re == 0.0 && z.im == 0.0 && exp < 0 {
                    Err(generic_error("Cannot raise zero complex number to negative power", None))
                } else {
                    Ok(Complex(z.powi(exp as i32)))
                }
            },

            Fraction(r) => {
                if r.is_zero() && exp < 0 {
                    Err(generic_error("Cannot raise zero fraction to negative power", None))
                } else if exp >= 0 {
                    Ok(Fraction(r.pow(exp as i32)))
                } else {
                    Ok(Fraction(r.recip().pow((-exp) as i32)))
                }
            },
        }
    }
}

// From<Number> implementation so that we can use into() method for Number
impl From<Number> for Arc<dyn Expr> {
    #[inline]
    fn from(num: Number) -> Self {
        intern_expr(Arc::new(num))
    }
}

impl From<&Number> for Arc<dyn Expr> {
    #[inline]
    fn from(num: &Number) -> Self {
        intern_expr(Arc::new(num.clone()))
    }
}

impl ExprInternal for Number {
    impl_expr_internal_methods!(Number, false);

    #[inline]
    fn hash_key(&self) -> String {
        match self {
            Number::Integer(n) => format!("Integer({})", n),
            Number::Real(f) => format!("Real({})", f),
            Number::Complex(z) => format!("Complex({})", z),
            Number::Fraction(r) => format!("Fraction({}/{})", r.numer(), r.denom()),
        }
    }

    #[inline]
    fn replace_expr_children(
        &self,
        _map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        _include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(self.clone_expr())
    }

    #[inline]
    fn is_exact_zero(&self) -> bool {
        self.approx_eq_number(&Number::Integer(0), Some(NumberTolerance::new(0.0, 0.0)))
    }
}

#[typetag::serde]
impl Expr for Number {
    impl_expr_common_methods!(true);

    #[allow(unused_variables)]
    #[inline]
    fn differentiate(
        &self,
        _s: &Arc<crate::perturbations::Perturbation>,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        Ok(Number::zero())
    }

    #[inline]
    fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_self_any(set, false) {
            Ok(Number::zero())
        } else {
            Ok(self.clone_expr())
        }
    }

    #[inline]
    fn retain(
        &self,
        set: &HashSet<Arc<dyn Expr>>,
        include_derivatives: bool,
    ) -> Result<Arc<dyn Expr>, TinnedError> {
        if self.match_self_any(set, include_derivatives) {
            Ok(self.clone_expr())
        } else {
            Ok(Number::zero())
        }
    }
}

impl PartialEq for Number {
    fn eq(&self, other: &Self) -> bool {
        self.approx_eq_number(other, None)
    }
}

impl Eq for Number {}

impl std::fmt::Display for Number {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Number::Integer(n) => write!(f, "{}", n),
            Number::Real(r) => write!(f, "{}", r),
            Number::Complex(z) => write!(f, "{} + {}i", z.re, z.im),
            Number::Fraction(r) => write!(f, "{}/{}", r.numer(), r.denom()),
        }
    }
}

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use rand::random_range;

    // Returns the `f64` that is `n` ULPs greater than the given `base` value.
    // If `n == 0`, returns the base itself.
    #[inline]
    pub fn ulps_up_f64(base: f64, n: u64) -> f64 {
        f64::from_bits(base.to_bits().wrapping_add(n))
    }

    // Returns the `f64` that is `n` ULPs less than the given `base` value.
    // Be careful near 0 — subnormals and negatives can behave differently.
    #[inline]
    pub fn ulps_down_f64(base: f64, n: u64) -> f64 {
        f64::from_bits(base.to_bits().wrapping_sub(n))
    }

    #[inline]
    pub fn make_number_i64(val_range: u32) -> Arc<dyn Expr> {
        let val: i64 = random_range(-(val_range as i64)..=val_range as i64);
        Number::from_i64(val)
    }

    #[inline]
    pub fn make_number_f64(val_range: u32) -> Arc<dyn Expr> {
        let val: f64 = random_range(-(val_range as f64)..=val_range as f64);
        Number::from_f64(val)
    }

    #[inline]
    pub fn make_number_complex(val_range: u32) -> Arc<dyn Expr> {
        let real: f64 = random_range(-(val_range as f64)..=val_range as f64);
        let imaginary: f64 = random_range(-(val_range as f64)..=val_range as f64);
        Number::from_complex(Complex64::new(real, imaginary))
    }

    #[inline]
    pub fn make_number_rational(val_range: u32) -> Arc<dyn Expr> {
        let numerator: i64 = random_range(-(val_range as i64)..=val_range as i64);
        let denominator: i64 = random_range(1..=val_range.max(1) as i64);
        Number::from_rational(Rational64::new(numerator, denominator))
    }
}

#[cfg(test)]
mod tests {
    use super::test_utils::*;
    use super::*;
    use crate::perturbations::perturbation::test_utils::make_perturbation_symbol;
    use crate::public::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};
    use num_integer::Integer;
    use rand::random_range;

    test_struct_safety!(Number);

    test_thread_interning!({
        let expr: Arc<dyn Expr> = Number::Integer(123).into();
        expr
    });

    // Basic structure and methods
    #[test]
    fn test_struct() {
        assert!(Number::Integer(0).is_zero(None));
        assert!(Number::Real(0.0).is_zero(None));
        assert!(Number::Complex(Complex64::new(0.0, 0.0)).is_zero(None));
        assert!(Number::Fraction(Rational64::new(0, 1)).is_zero(None));
        assert!(!Number::Integer(1).is_zero(None));
        assert!(!Number::Real(1.0).is_zero(None));
        assert!(!Number::Complex(Complex64::new(1.0, 0.0)).is_zero(None));
        assert!(!Number::Fraction(Rational64::new(1, 1)).is_zero(None));

        assert!(Number::Integer(1).is_one(None));
        assert!(Number::Real(1.0).is_one(None));
        assert!(Number::Complex(Complex64::new(1.0, 0.0)).is_one(None));
        assert!(Number::Fraction(Rational64::new(1, 1)).is_one(None));
        assert!(!Number::Integer(0).is_one(None));
        assert!(!Number::Real(0.0).is_one(None));
        assert!(!Number::Complex(Complex64::new(0.0, 0.0)).is_one(None));
        assert!(!Number::Fraction(Rational64::new(0, 1)).is_one(None));

        let mut real = Number::Real(-0.0);
        let mut cmplx = Number::Complex(Complex64::new(0.0, -0.0));

        assert!(real.is_zero(None));
        assert!(cmplx.is_zero(None));

        let abs_error: f64 = 0.0;
        let rel_error: f64 = 1e-12;
        let tol = NumberTolerance::new(abs_error, rel_error);
        let tight_tol = NumberTolerance::new(abs_error, 0.25 * rel_error);
        real = Number::Real(1.0 + 0.5 * rel_error);
        cmplx = Number::Complex(Complex64::new(1.0 + 0.5 * rel_error, -0.0));

        assert!(!real.is_one(Some(tight_tol.clone())));
        assert!(!cmplx.is_one(Some(tight_tol.clone())));
        assert!(real.is_one(Some(tol.clone())));
        assert!(cmplx.is_one(Some(tol.clone())));

        real = Number::Real(1.0 + 5.0 * rel_error);
        cmplx = Number::Complex(Complex64::new(1.0 + 5.0 * rel_error, -0.0));

        assert!(!real.is_one(Some(tol.clone())));
        assert!(!cmplx.is_one(Some(tol.clone())));

        let n1: i64 = random_range(-100..=100);
        let n2: i64 = random_range(1..=100);
        let f1: f64 = random_range(-100.0..=100.0);
        let f2: f64 = random_range(-100.0..=100.0);

        assert_eq!(Number::Integer(n1), Number::Integer(n1));
        assert_eq!(Number::Integer(n1), Number::Real(n1 as f64));
        assert_eq!(Number::Real(f1), Number::Real(ulps_up_f64(f1, 3)));
        assert_eq!(Number::Real(f1), Number::Real(ulps_down_f64(f1, 3)));
        assert_eq!(
            Number::Complex(Complex64::new(f1, f2)),
            Number::Complex(Complex64::new(ulps_up_f64(f1, 3), ulps_down_f64(f2, 3)))
        );
        assert_eq!(
            Number::Fraction(Rational64::new(n1, n2)),
            Number::Fraction(Rational64::new(n1, n2))
        );
    }

    #[test]
    fn test_negate() {
        let n1: i64 = random_range(-100..=100);
        assert_eq!(Number::Integer(n1).negate(), Number::Integer(-n1));

        let f1: f64 = random_range(-100.0..=100.0);
        assert_eq!(Number::Real(f1).negate(), Number::Real(-f1));

        let f2: f64 = random_range(-100.0..=100.0);
        assert_eq!(
            Number::Complex(Complex64::new(f1, f2)).negate(),
            Number::Complex(Complex64::new(-f1, -f2))
        );

        let n2: i64 = random_range(1..=100);
        assert_eq!(
            Number::Fraction(Rational64::new(n1, n2)).negate(),
            Number::Fraction(Rational64::new(-n1, n2))
        );
    }

    #[test]
    fn test_conjugate() {
        let n1: i64 = random_range(-100..=100);
        assert_eq!(Number::Integer(n1).conjugate(), Number::Integer(n1));

        let f1: f64 = random_range(-100.0..=100.0);
        assert_eq!(Number::Real(f1).conjugate(), Number::Real(f1));

        let f2: f64 = random_range(-100.0..=100.0);
        assert_eq!(
            Number::Complex(Complex64::new(f1, f2)).conjugate(),
            Number::Complex(Complex64::new(f1, -f2))
        );

        let n2: i64 = random_range(1..=100);
        assert_eq!(
            Number::Fraction(Rational64::new(n1, n2)).conjugate(),
            Number::Fraction(Rational64::new(n1, n2))
        );
    }

    #[test]
    fn test_add() {
        let n1: i64 = random_range(-100..=100);
        let n2: i64 = random_range(1..=100);
        let f1: f64 = random_range(-100.0..=100.0);
        let f2: f64 = random_range(-100.0..=100.0);

        assert_eq!(Number::Integer(n1).add(&Number::Integer(n2)), Number::Integer(n1 + n2));
        assert_eq!(Number::Integer(n1).add(&Number::Real(f1)), Number::Real(n1 as f64 + f1));
        assert_eq!(
            Number::Integer(n1).add(&Number::Complex(Complex64::new(f1, f2))),
            Number::Complex(Complex64::new(n1 as f64 + f1, f2))
        );
        assert_eq!(
            Number::Integer(n1).add(&Number::Fraction(Rational64::new(n1, n2))),
            Number::Fraction(Rational64::new(n1 * (n2 + 1), n2))
        );

        assert_eq!(Number::Real(f1).add(&Number::Integer(n1)), Number::Real(f1 + n1 as f64));
        assert_eq!(Number::Real(f1).add(&Number::Real(f2)), Number::Real(f1 + f2));
        assert_eq!(
            Number::Real(f1).add(&Number::Complex(Complex64::new(f1, f2))),
            Number::Complex(Complex64::new(f1 + f1, f2))
        );
        assert_eq!(
            Number::Real(f1).add(&Number::Fraction(Rational64::new(n1, n2))),
            Number::Real(f1 + (n1 as f64) / (n2 as f64))
        );

        assert_eq!(
            Number::Complex(Complex64::new(f1, f2)).add(&Number::Integer(n1)),
            Number::Complex(Complex64::new(f1 + n1 as f64, f2))
        );
        assert_eq!(
            Number::Complex(Complex64::new(f1, f2)).add(&Number::Real(f1)),
            Number::Complex(Complex64::new(f1 + f1, f2))
        );
        assert_eq!(
            Number::Complex(Complex64::new(f1, f2)).add(&Number::Complex(Complex64::new(f1, f2))),
            Number::Complex(Complex64::new(f1 + f1, f2 + f2))
        );
        assert_eq!(
            Number::Complex(Complex64::new(f1, f2)).add(&Number::Fraction(Rational64::new(n1, n2))),
            Number::Complex(Complex64::new(f1 + (n1 as f64) / (n2 as f64), f2))
        );

        assert_eq!(
            Number::Fraction(Rational64::new(n1, n2)).add(&Number::Integer(n1)),
            Number::Fraction(Rational64::new(n1 * (n2 + 1), n2))
        );
        assert_eq!(
            Number::Fraction(Rational64::new(n1, n2)).add(&Number::Real(f1)),
            Number::Real((n1 as f64) / (n2 as f64) + f1)
        );
        assert_eq!(
            Number::Fraction(Rational64::new(n1, n2)).add(&Number::Complex(Complex64::new(f1, f2))),
            Number::Complex(Complex64::new((n1 as f64) / (n2 as f64) + f1, f2))
        );
        assert_eq!(
            Number::Fraction(Rational64::new(n1, n2))
                .add(&Number::Fraction(Rational64::new(n1, n2))),
            Number::Fraction(Rational64::new(n1 + n1, n2))
        );
    }

    #[test]
    fn test_mul() {
        let n1: i64 = random_range(-100..=100);
        let n2: i64 = random_range(1..=100);
        let f1: f64 = random_range(-100.0..=100.0);
        let f2: f64 = random_range(-100.0..=100.0);

        assert_eq!(Number::Integer(n1).mul(&Number::Integer(n2)), Number::Integer(n1 * n2));
        assert_eq!(Number::Integer(n1).mul(&Number::Real(f1)), Number::Real(n1 as f64 * f1));
        assert_eq!(
            Number::Integer(n1).mul(&Number::Complex(Complex64::new(f1, f2))),
            Number::Complex(Complex64::new(n1 as f64 * f1, n1 as f64 * f2))
        );
        assert_eq!(
            Number::Integer(n1).mul(&Number::Fraction(Rational64::new(n1, n2))),
            Number::Fraction(Rational64::new(n1 * n1, n2))
        );

        assert_eq!(Number::Real(f1).mul(&Number::Integer(n1)), Number::Real(f1 * n1 as f64));
        assert_eq!(Number::Real(f1).mul(&Number::Real(f2)), Number::Real(f1 * f2));
        assert_eq!(
            Number::Real(f1).mul(&Number::Complex(Complex64::new(f1, f2))),
            Number::Complex(Complex64::new(f1 * f1, f1 * f2))
        );
        assert_eq!(
            Number::Real(f1).mul(&Number::Fraction(Rational64::new(n1, n2))),
            Number::Real(f1 * (n1 as f64) / (n2 as f64))
        );

        assert_eq!(
            Number::Complex(Complex64::new(f1, f2)).mul(&Number::Integer(n1)),
            Number::Complex(Complex64::new(f1 * n1 as f64, f2 * n1 as f64))
        );
        assert_eq!(
            Number::Complex(Complex64::new(f1, f2)).mul(&Number::Real(f1)),
            Number::Complex(Complex64::new(f1 * f1, f2 * f1))
        );
        assert_eq!(
            Number::Complex(Complex64::new(f1, f2)).mul(&Number::Complex(Complex64::new(f1, f2))),
            Number::Complex(Complex64::new(f1 * f1 - f2 * f2, f1 * f2 + f2 * f1))
        );
        assert_eq!(
            Number::Complex(Complex64::new(f1, f2)).mul(&Number::Fraction(Rational64::new(n1, n2))),
            Number::Complex(Complex64::new(
                f1 * (n1 as f64) / (n2 as f64),
                f2 * (n1 as f64) / (n2 as f64)
            ))
        );

        assert_eq!(
            Number::Fraction(Rational64::new(n1, n2)).mul(&Number::Integer(n1)),
            Number::Fraction(Rational64::new(n1 * n1, n2))
        );
        assert_eq!(
            Number::Fraction(Rational64::new(n1, n2)).mul(&Number::Real(f1)),
            Number::Real((n1 as f64) / (n2 as f64) * f1)
        );
        assert_eq!(
            Number::Fraction(Rational64::new(n1, n2)).mul(&Number::Complex(Complex64::new(f1, f2))),
            Number::Complex(Complex64::new(
                (n1 as f64) / (n2 as f64) * f1,
                (n1 as f64) / (n2 as f64) * f2
            ))
        );
        assert_eq!(
            Number::Fraction(Rational64::new(n1, n2))
                .mul(&Number::Fraction(Rational64::new(n1, n2))),
            Number::Fraction(Rational64::new(n1 * n1, n2 * n2))
        );
    }

    #[test]
    fn test_pow_i64() {
        assert_eq!(Number::Integer(2).pow_i64(3).unwrap(), Number::Integer(8));
        assert_eq!(
            Number::Integer(2).pow_i64(-2).unwrap(),
            Number::Fraction(Rational64::new(1, 4))
        );
        assert_eq!(Number::Real(2.0).pow_i64(2).unwrap(), Number::Real(4.0));
        assert_eq!(
            Number::Complex(Complex64::new(0.0, 1.0)).pow_i64(2).unwrap(),
            Number::Complex(Complex64::new(-1.0, 0.0))
        );
        assert_eq!(
            Number::Complex(Complex64::new(0.0, 1.0)).pow_i64(-3).unwrap(),
            Number::Complex(Complex64::new(0.0, 1.0))
        );
        assert_eq!(
            Number::Fraction(Rational64::new(1, 2)).pow_i64(2).unwrap(),
            Number::Fraction(Rational64::new(1, 4))
        );
    }

    // Implementation for ExprInternal and Expr
    #[test]
    fn test_impl_expr() {
        let n1: i64 = random_range(-100..=100);
        let n2: i64 = random_range(2..=100);
        let f1: f64 = random_range(-100.0..=100.0);
        let f2: f64 = random_range(-100.0..=100.0);

        let int: Arc<dyn Expr> = Number::Integer(n1).into();
        let real: Arc<dyn Expr> = Number::Real(f1).into();
        let cmplx: Arc<dyn Expr> = Number::Complex(Complex64::new(f1, f2)).into();
        let frac: Arc<dyn Expr> = Number::Fraction(Rational64::new(n1, n2)).into();

        assert_eq!(&int, &Number::from_i64(n1));
        assert_eq!(&real, &Number::from_f64(f1));
        assert_eq!(&cmplx, &Number::from_complex(Complex64::new(f1, f2)));
        assert_eq!(&frac, &Number::from_rational(Rational64::new(n1, n2)));

        assert_eq!(format!("{}", int), n1.to_string());
        assert_eq!(format!("{}", real), f1.to_string());
        assert_eq!(format!("{}", cmplx), format!("{f1} + {f2}i"));

        let g = n1.gcd(&n2);
        let new_n1 = n1 / g;
        let new_n2 = n2 / g;

        assert_eq!(format!("{}", frac), format!("{new_n1}/{new_n2}"));

        assert_eq!(int.hash_key(), format!("Integer({n1})"));
        assert_eq!(real.hash_key(), format!("Real({f1})"));

        if f2 >= 0.0 {
            assert_eq!(cmplx.hash_key(), format!("Complex({f1}+{f2}i)"));
        } else {
            assert_eq!(cmplx.hash_key(), format!("Complex({f1}{f2}i)"));
        }

        assert_eq!(frac.hash_key(), format!("Fraction({new_n1}/{new_n2})"));

        assert!(int.total_order() == 0);
        assert!(real.total_order() == 0);
        assert!(cmplx.total_order() == 0);
        assert!(frac.total_order() == 0);

        assert!(int.deep_eq_superchains(&int));
        assert!(!int.deep_eq_superchains(&real));
        assert!(!int.deep_eq_superchains(&cmplx));
        assert!(!int.deep_eq_superchains(&frac));
        assert!(!real.deep_eq_superchains(&int));
        assert!(real.deep_eq_superchains(&real));
        assert!(!real.deep_eq_superchains(&cmplx));
        assert!(!real.deep_eq_superchains(&frac));
        assert!(!cmplx.deep_eq_superchains(&int));
        assert!(!cmplx.deep_eq_superchains(&real));
        assert!(cmplx.deep_eq_superchains(&cmplx));
        assert!(!cmplx.deep_eq_superchains(&frac));
        assert!(!frac.deep_eq_superchains(&int));
        assert!(!frac.deep_eq_superchains(&real));
        assert!(!frac.deep_eq_superchains(&cmplx));
        assert!(frac.deep_eq_superchains(&frac));

        assert!(int.eq_by_superchains(&int));
        assert!(!int.eq_by_superchains(&real));
        assert!(!int.eq_by_superchains(&cmplx));
        assert!(!int.eq_by_superchains(&frac));
        assert!(!real.eq_by_superchains(&int));
        assert!(real.eq_by_superchains(&real));
        assert!(!real.eq_by_superchains(&cmplx));
        assert!(!real.eq_by_superchains(&frac));
        assert!(!cmplx.eq_by_superchains(&int));
        assert!(!cmplx.eq_by_superchains(&real));
        assert!(cmplx.eq_by_superchains(&cmplx));
        assert!(!cmplx.eq_by_superchains(&frac));
        assert!(!frac.eq_by_superchains(&int));
        assert!(!frac.eq_by_superchains(&real));
        assert!(!frac.eq_by_superchains(&cmplx));
        assert!(frac.eq_by_superchains(&frac));

        //FIXME: add the following tests
        //assert_eq!(&int.replace_expr_self().unwrap(), &int);
        //assert_eq!(&real.replace_expr_self().unwrap(), &real);
        //assert_eq!(&cmplx.replace_expr_self().unwrap(), &cmplx);
        //assert_eq!(&frac.replace_expr_self().unwrap(), &frac);

        //replace_expr_self
        //replace_expr_children
        //retain_expr_fields

        assert!(int.is_scalar());
        assert!(real.is_scalar());
        assert!(cmplx.is_scalar());
        assert!(frac.is_scalar());
    }

    #[test]
    fn test_differentiation() {
        let int = make_number_i64(100u32);
        let real = make_number_f64(100u32);
        let cmplx = make_number_complex(100u32);
        let frac = make_number_rational(100u32);

        let p = make_perturbation_symbol(4u32, 4u32);

        assert_eq!(&int.differentiate(&p).unwrap(), &Number::zero());
        assert_eq!(&real.differentiate(&p).unwrap(), &Number::zero());
        assert_eq!(&cmplx.differentiate(&p).unwrap(), &Number::zero());
        assert_eq!(&frac.differentiate(&p).unwrap(), &Number::zero());
    }

    // Test serialization and deserialization via `serde_json`
    #[test]
    fn test_serialization() {
        let mut n = make_number_i64(100u32);
        let mut json = serde_json::to_string(&n).unwrap();
        let mut recovered: Arc<dyn Expr> = serde_json::from_str(&json).unwrap();
        assert_eq!(&n, &recovered);

        n = make_number_f64(100u32);
        json = serde_json::to_string(&n).unwrap();
        recovered = serde_json::from_str(&json).unwrap();
        assert_eq!(&n, &recovered);

        n = make_number_complex(100u32);
        json = serde_json::to_string(&n).unwrap();
        recovered = serde_json::from_str(&json).unwrap();
        assert_eq!(&n, &recovered);

        n = make_number_rational(100u32);
        json = serde_json::to_string(&n).unwrap();
        recovered = serde_json::from_str(&json).unwrap();
        assert_eq!(&n, &recovered);
    }

    // Test internal and public utils
    #[test]
    fn test_utils() {
        let n1: i64 = random_range(-100..=100);
        let n2: i64 = random_range(1..=100);
        let f1: f64 = random_range(-100.0..=100.0);
        let f2: f64 = random_range(-100.0..=100.0);

        let int = Number::from_i64(n1);
        let real = Number::from_f64(f1);
        let cmplx = Number::from_complex(Complex64::new(f1, f2));
        let frac = Number::from_rational(Rational64::new(n1, n2));

        assert!(is_expr_type::<Number>(&int));
        assert!(is_expr_type::<Number>(&real));
        assert!(is_expr_type::<Number>(&cmplx));
        assert!(is_expr_type::<Number>(&frac));

        assert_eq!(is_zero_expr(&int, None), n1 == 0);
        assert_eq!(is_zero_expr(&real, None), f1 == 0.0);
        assert_eq!(is_zero_expr(&cmplx, None), f1 == 0.0 && f2 == 0.0);
        assert_eq!(is_zero_expr(&frac, None), n1 == 0);

        assert_eq!(is_one_expr(&int, None), n1 == 1);
        assert_eq!(is_one_expr(&real, None), f1 == 1.0);
        assert_eq!(is_one_expr(&cmplx, None), f1 == 1.0 && f2 == 0.0);
        assert_eq!(is_one_expr(&frac, None), n1 == n2);

        let mut num = downcast_from_arc::<Number>(&int).unwrap();
        assert_eq!(num, &Number::Integer(n1));

        num = downcast_from_arc::<Number>(&real).unwrap();
        assert_eq!(num, &Number::Real(f1));

        num = downcast_from_arc::<Number>(&cmplx).unwrap();
        assert_eq!(num, &Number::Complex(Complex64::new(f1, f2)));

        num = downcast_from_arc::<Number>(&frac).unwrap();
        assert_eq!(num, &Number::Fraction(Rational64::new(n1, n2)));

        let int1 = Number::Integer(n1).into();
        let int2 = Number::from_i64(n1 + 1);

        assert!(Arc::ptr_eq(&int, &int1));
        assert!(!Arc::ptr_eq(&int, &int2));

        let real1 = Number::Real(f1).into();
        let real2 = Number::from_f64(f1 + 1.0);

        assert!(Arc::ptr_eq(&real, &real1));
        assert!(!Arc::ptr_eq(&real, &real2));

        let cmplx1 = Number::Complex(Complex64::new(f1, f2)).into();
        let cmplx2 = Number::from_complex(Complex64::new(f1 + 1.0, f2));

        assert!(Arc::ptr_eq(&cmplx, &cmplx1));
        assert!(!Arc::ptr_eq(&cmplx, &cmplx2));

        let frac1 = Number::Fraction(Rational64::new(n1, n2)).into();
        let frac2 = Number::from_rational(Rational64::new(n1 + 1, n2));

        assert!(Arc::ptr_eq(&frac, &frac1));
        assert!(!Arc::ptr_eq(&frac, &frac2));
    }
}
