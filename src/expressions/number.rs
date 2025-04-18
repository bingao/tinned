use std::sync::Arc;

use float_cmp::approx_eq;
use num_complex::Complex64;
use num_rational::Rational64;
use num_traits::ToPrimitive;

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

            (Real(a), Fraction(b)) => {
                Real(*a + b.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", b)))
            },
            (Fraction(a), Real(b)) => {
                Real(a.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", a)) + *b)
            },

            (Real(a), Complex(b)) => Complex(Complex64::new(*a, 0.0) + *b),
            (Complex(a), Real(b)) => Complex(*a + Complex64::new(*b, 0.0)),

            (Fraction(a), Complex(b)) => Complex(
                Complex64::new(
                    a.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", a)),
                    0.0,
                ) + *b,
            ),
            (Complex(a), Fraction(b)) => Complex(
                *a + Complex64::new(
                    b.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", b)),
                    0.0,
                ),
            ),

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

            (Real(a), Fraction(b)) => {
                Real(*a * b.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", b)))
            },
            (Fraction(a), Real(b)) => {
                Real(a.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", a)) * *b)
            },

            (Real(a), Complex(b)) => Complex(Complex64::new(*a, 0.0) * *b),
            (Complex(a), Real(b)) => Complex(*a * Complex64::new(*b, 0.0)),

            (Fraction(a), Complex(b)) => Complex(
                Complex64::new(
                    a.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", a)),
                    0.0,
                ) * *b,
            ),
            (Complex(a), Fraction(b)) => Complex(
                *a * Complex64::new(
                    b.to_f64().unwrap_or_else(|| panic!("Failed to convert {} to f64", b)),
                    0.0,
                ),
            ),

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

impl From<&Number> for Arc<dyn Expr> {
    #[inline]
    fn from(num: &Number) -> Self {
        crate::utils::intern_expr(Arc::new(num.clone()))
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
            (Number::Real(a), Number::Real(b)) => {
                approx_eq!(f64, *a, *b, ulps = 4)
            },
            (Number::Complex(a), Number::Complex(b)) => {
                approx_eq!(f64, a.re, b.re, ulps = 4) && approx_eq!(f64, a.im, b.im, ulps = 4)
            },
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
    use crate::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};
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
        assert!(Number::Integer(0).is_zero());
        assert!(Number::Real(0.0).is_zero());
        assert!(Number::Complex(Complex64::new(0.0, 0.0)).is_zero());
        assert!(Number::Fraction(Rational64::new(0, 1)).is_zero());
        assert!(!Number::Integer(1).is_zero());
        assert!(!Number::Real(1.0).is_zero());
        assert!(!Number::Complex(Complex64::new(1.0, 0.0)).is_zero());
        assert!(!Number::Fraction(Rational64::new(1, 1)).is_zero());

        assert!(Number::Integer(1).is_one());
        assert!(Number::Real(1.0).is_one());
        assert!(Number::Complex(Complex64::new(1.0, 0.0)).is_one());
        assert!(Number::Fraction(Rational64::new(1, 1)).is_one());
        assert!(!Number::Integer(0).is_one());
        assert!(!Number::Real(0.0).is_one());
        assert!(!Number::Complex(Complex64::new(0.0, 0.0)).is_one());
        assert!(!Number::Fraction(Rational64::new(0, 1)).is_one());

        let n1: i64 = random_range(-100..=100);
        let f1: f64 = random_range(-100.0..=100.0);

        assert_eq!(Number::Integer(n1).hash_key(), format!("Integer({n1})"));
        assert_eq!(Number::Real(f1).hash_key(), format!("Real({f1})"));

        let f2: f64 = random_range(-100.0..=100.0);
        if f2 >= 0.0 {
            assert_eq!(
                Number::Complex(Complex64::new(f1, f2)).hash_key(),
                format!("Complex({f1}+{f2}i)")
            );
        } else {
            assert_eq!(
                Number::Complex(Complex64::new(f1, f2)).hash_key(),
                format!("Complex({f1}{f2}i)")
            );
        }

        let n2: i64 = random_range(1..=100);
        let g = n1.gcd(&n2);
        let new_n1 = n1 / g;
        let new_n2 = n2 / g;
        assert_eq!(
            Number::Fraction(Rational64::new(n1, n2)).hash_key(),
            format!("Fraction({new_n1}/{new_n2})")
        );

        assert!(Number::Integer(n1).is_scalar());
        assert!(Number::Real(f1).is_scalar());
        assert!(Number::Complex(Complex64::new(f1, f2)).is_scalar());
        assert!(Number::Fraction(Rational64::new(n1, n2)).is_scalar());

        assert_eq!(Number::Integer(n1), Number::Integer(n1));
        assert_ne!(Number::Integer(n1), Number::Real(n1 as f64)); // they are different variants
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

        assert_eq!(format!("{}", Number::Integer(n1)), n1.to_string());
        assert_eq!(format!("{}", Number::Real(f1)), f1.to_string());
        assert_eq!(format!("{}", Number::Complex(Complex64::new(f1, f2))), format!("{f1} + {f2}i"));
        assert_eq!(
            format!("{}", Number::Fraction(Rational64::new(n1, n2))),
            format!("{new_n1}/{new_n2}")
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

    // Implementation for Expr
    #[test]
    fn test_impl_expr() {
        let n1: i64 = random_range(-100..=100);
        let n2: i64 = random_range(1..=100);
        let f1: f64 = random_range(-100.0..=100.0);
        let f2: f64 = random_range(-100.0..=100.0);

        let int: Arc<dyn Expr> = Number::Integer(n1).into();
        let real: Arc<dyn Expr> = Number::Real(f1).into();
        let cmplx: Arc<dyn Expr> = Number::Complex(Complex64::new(f1, f2)).into();
        let frac: Arc<dyn Expr> = Number::Fraction(Rational64::new(n1, n2)).into();

        assert_eq!(int.hash_key(), format!("Integer({n1})"));
        assert_eq!(real.hash_key(), format!("Real({f1})"));

        if f2 >= 0.0 {
            assert_eq!(cmplx.hash_key(), format!("Complex({f1}+{f2}i)"));
        } else {
            assert_eq!(cmplx.hash_key(), format!("Complex({f1}{f2}i)"));
        }

        let g = n1.gcd(&n2);
        let new_n1 = n1 / g;
        let new_n2 = n2 / g;
        assert_eq!(frac.hash_key(), format!("Fraction({new_n1}/{new_n2})"));

        assert!(int.is_scalar());
        assert!(real.is_scalar());
        assert!(cmplx.is_scalar());
        assert!(frac.is_scalar());

        assert_eq!(&int, &Number::from_i64(n1));
        assert_eq!(&real, &Number::from_f64(f1));
        assert_eq!(&cmplx, &Number::from_complex(Complex64::new(f1, f2)));
        assert_eq!(&frac, &Number::from_rational(Rational64::new(n1, n2)));

        assert_ne!(&int, &real.clone());
        assert_ne!(&int, &cmplx.clone());
        assert_ne!(&int, &frac.clone());
        assert_ne!(&real, &cmplx.clone());
        assert_ne!(&real, &frac.clone());
        assert_ne!(&cmplx, &frac.clone());

        assert_eq!(format!("{}", int), n1.to_string());
        assert_eq!(format!("{}", real), f1.to_string());
        assert_eq!(format!("{}", cmplx), format!("{f1} + {f2}i"));
        assert_eq!(format!("{}", frac), format!("{new_n1}/{new_n2}"));
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

    // Test utils: interning, downcast, type and identity check
    #[test]
    fn test_utils() {
        let n1: i64 = random_range(-100..=100);
        let int = Number::from_i64(n1);
        let int1 = Number::Integer(n1).into();
        let int2 = Number::from_i64(n1 + 1);

        assert!(Arc::ptr_eq(&int, &int1));
        assert!(!Arc::ptr_eq(&int, &int2));

        let f1: f64 = random_range(-100.0..=100.0);
        let real = Number::from_f64(f1);
        let real1 = Number::Real(f1).into();
        let real2 = Number::from_f64(f1 + 1.0);

        assert!(Arc::ptr_eq(&real, &real1));
        assert!(!Arc::ptr_eq(&real, &real2));

        let f2: f64 = random_range(-100.0..=100.0);
        let cmplx = Number::from_complex(Complex64::new(f1, f2));
        let cmplx1 = Number::Complex(Complex64::new(f1, f2)).into();
        let cmplx2 = Number::from_complex(Complex64::new(f1 + 1.0, f2));

        assert!(Arc::ptr_eq(&cmplx, &cmplx1));
        assert!(!Arc::ptr_eq(&cmplx, &cmplx2));

        let n2: i64 = random_range(1..=100);
        let frac = Number::from_rational(Rational64::new(n1, n2));
        let frac1 = Number::Fraction(Rational64::new(n1, n2)).into();
        let frac2 = Number::from_rational(Rational64::new(n1 + 1, n2));

        assert!(Arc::ptr_eq(&frac, &frac1));
        assert!(!Arc::ptr_eq(&frac, &frac2));

        let mut num = downcast_from_arc::<Number>(&int).unwrap();
        assert_eq!(num, &Number::Integer(n1));

        num = downcast_from_arc::<Number>(&real).unwrap();
        assert_eq!(num, &Number::Real(f1));

        num = downcast_from_arc::<Number>(&cmplx).unwrap();
        assert_eq!(num, &Number::Complex(Complex64::new(f1, f2)));

        num = downcast_from_arc::<Number>(&frac).unwrap();
        assert_eq!(num, &Number::Fraction(Rational64::new(n1, n2)));

        assert!(is_expr_type::<Number>(&int));
        assert!(is_expr_type::<Number>(&real));
        assert!(is_expr_type::<Number>(&cmplx));
        assert!(is_expr_type::<Number>(&frac));

        assert_eq!(is_zero_expr(&int), n1 == 0);
        assert_eq!(is_zero_expr(&real), f1 == 0.0);
        assert_eq!(is_zero_expr(&cmplx), f1 == 0.0 && f2 == 0.0);
        assert_eq!(is_zero_expr(&frac), n1 == 0);

        assert_eq!(is_one_expr(&int), n1 == 1);
        assert_eq!(is_one_expr(&real), f1 == 1.0);
        assert_eq!(is_one_expr(&cmplx), f1 == 1.0 && f2 == 0.0);
        assert_eq!(is_one_expr(&frac), n1 == n2);
    }
}
