use std::sync::Arc;

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

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use num_complex::Complex64;
    use num_rational::Rational64;

    use crate::core::Expr;
    use crate::expressions::Number;
    use crate::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};

    // Basic structure and methods
    #[test]
    fn test_struct() {
        assert!(Number::Integer(0).is_zero());
        assert!(Number::Real(0.0).is_zero());
        assert!(Number::Fraction(Rational64::new(0, 1)).is_zero());
        assert!(Number::Complex(Complex64::new(0.0, 0.0)).is_zero());
        assert!(!Number::Integer(1).is_zero());

        assert!(Number::Integer(1).is_one());
        assert!(Number::Real(1.0).is_one());
        assert!(Number::Fraction(Rational64::new(1, 1)).is_one());
        assert!(Number::Complex(Complex64::new(1.0, 0.0)).is_one());
        assert!(!Number::Integer(0).is_one());

        assert_eq!(Number::Integer(1).hash_key(), "Integer(1)");
        assert_eq!(Number::Real(0.0).hash_key(), "Real(0)");
        assert_eq!(Number::Complex(Complex64::new(1.0, -1.0)).hash_key(), "Complex(1-1i)");
        assert_eq!(Number::Fraction(Rational64::new(3, 4)).hash_key(), "Fraction(3/4)");

        assert!(Number::Integer(1).is_scalar());
        assert!(Number::Real(0.0).is_scalar());
        assert!(Number::Complex(Complex64::new(1.0, -1.0)).is_scalar());
        assert!(Number::Fraction(Rational64::new(3, 4)).is_scalar());

        assert_eq!(Number::Integer(5), Number::Integer(5));
        assert_ne!(Number::Integer(5), Number::Real(5.0)); // they are different variants
        assert_eq!(Number::Real(5.0), Number::Real(5.0));
        assert_eq!(
            Number::Complex(Complex64::new(2.0, 3.0)),
            Number::Complex(Complex64::new(2.0, 3.0))
        );
        assert_eq!(
            Number::Fraction(Rational64::new(1, 2)),
            Number::Fraction(Rational64::new(1, 2))
        );

        assert_eq!(format!("{}", Number::Integer(42)), "42");
        assert_eq!(format!("{}", Number::Real(3.14)), "3.14");
        assert_eq!(format!("{}", Number::Fraction(Rational64::new(3, 2))), "3/2");
        assert_eq!(format!("{}", Number::Complex(Complex64::new(1.0, 2.0))), "1 + 2i");
    }

    #[test]
    fn test_conjugate() {
        let c = Number::Complex(Complex64::new(3.0, 4.0));
        if let Number::Complex(z) = c.conjugate() {
            assert_eq!(z.re, 3.0);
            assert_eq!(z.im, -4.0);
        } else {
            panic!("Conjugate of complex should remain complex");
        }

        assert_eq!(Number::Integer(5).conjugate(), Number::Integer(5));
    }

    #[test]
    fn test_add() {
        assert_eq!(Number::Integer(2).add(&Number::Integer(3)), Number::Integer(5));
        assert_eq!(Number::Integer(2).add(&Number::Real(3.0)), Number::Real(5.0));
        assert_eq!(
            Number::Integer(2).add(&Number::Complex(Complex64::new(1.0, -3.0))),
            Number::Complex(Complex64::new(3.0, -3.0))
        );
        assert_eq!(
            Number::Integer(2).add(&Number::Fraction(Rational64::new(3, 4))),
            Number::Fraction(Rational64::new(11, 4))
        );

        assert_eq!(Number::Real(2.5).add(&Number::Integer(3)), Number::Real(5.5));
        assert_eq!(Number::Real(2.5).add(&Number::Real(3.0)), Number::Real(5.5));
        assert_eq!(
            Number::Real(2.5).add(&Number::Complex(Complex64::new(1.0, -3.0))),
            Number::Complex(Complex64::new(3.5, -3.0))
        );
        assert_eq!(
            Number::Real(2.5).add(&Number::Fraction(Rational64::new(3, 4))),
            Number::Real(3.25)
        );

        assert_eq!(
            Number::Complex(Complex64::new(2.0, 3.0)).add(&Number::Integer(3)),
            Number::Complex(Complex64::new(5.0, 3.0))
        );
        assert_eq!(
            Number::Complex(Complex64::new(2.0, 3.0)).add(&Number::Real(3.0)),
            Number::Complex(Complex64::new(5.0, 3.0))
        );
        assert_eq!(
            Number::Complex(Complex64::new(2.0, 3.0))
                .add(&Number::Complex(Complex64::new(1.0, -3.0))),
            Number::Complex(Complex64::new(3.0, 0.0))
        );
        assert_eq!(
            Number::Complex(Complex64::new(2.0, 3.0)).add(&Number::Fraction(Rational64::new(3, 4))),
            Number::Complex(Complex64::new(2.75, 3.0))
        );

        assert_eq!(
            Number::Fraction(Rational64::new(1, 2)).add(&Number::Integer(3)),
            Number::Fraction(Rational64::new(7, 2))
        );
        assert_eq!(
            Number::Fraction(Rational64::new(1, 2)).add(&Number::Real(3.0)),
            Number::Real(3.5)
        );
        assert_eq!(
            Number::Fraction(Rational64::new(1, 2))
                .add(&Number::Complex(Complex64::new(1.0, -3.0))),
            Number::Complex(Complex64::new(1.5, -3.0))
        );
        assert_eq!(
            Number::Fraction(Rational64::new(1, 2)).add(&Number::Fraction(Rational64::new(3, 4))),
            Number::Fraction(Rational64::new(5, 4))
        );
    }

    #[test]
    fn test_mul() {
        assert_eq!(Number::Integer(2).mul(&Number::Integer(3)), Number::Integer(6));
        assert_eq!(Number::Integer(2).mul(&Number::Real(3.0)), Number::Real(6.0));
        assert_eq!(
            Number::Integer(2).mul(&Number::Complex(Complex64::new(1.0, -3.0))),
            Number::Complex(Complex64::new(2.0, -6.0))
        );
        assert_eq!(
            Number::Integer(2).mul(&Number::Fraction(Rational64::new(3, 4))),
            Number::Fraction(Rational64::new(3, 2))
        );

        assert_eq!(Number::Real(2.5).mul(&Number::Integer(3)), Number::Real(7.5));
        assert_eq!(Number::Real(2.5).mul(&Number::Real(3.0)), Number::Real(7.5));
        assert_eq!(
            Number::Real(2.5).mul(&Number::Complex(Complex64::new(1.0, -3.0))),
            Number::Complex(Complex64::new(2.5, -7.5))
        );
        assert_eq!(
            Number::Real(2.5).mul(&Number::Fraction(Rational64::new(3, 4))),
            Number::Real(1.875)
        );

        assert_eq!(
            Number::Complex(Complex64::new(2.0, 3.0)).mul(&Number::Integer(3)),
            Number::Complex(Complex64::new(6.0, 9.0))
        );
        assert_eq!(
            Number::Complex(Complex64::new(2.0, 3.0)).mul(&Number::Real(3.0)),
            Number::Complex(Complex64::new(6.0, 9.0))
        );
        assert_eq!(
            Number::Complex(Complex64::new(2.0, 3.0))
                .mul(&Number::Complex(Complex64::new(1.0, -3.0))),
            Number::Complex(Complex64::new(11.0, -3.0))
        );
        assert_eq!(
            Number::Complex(Complex64::new(2.0, 3.0)).mul(&Number::Fraction(Rational64::new(3, 4))),
            Number::Complex(Complex64::new(1.5, 2.25))
        );

        assert_eq!(
            Number::Fraction(Rational64::new(1, 2)).mul(&Number::Integer(3)),
            Number::Fraction(Rational64::new(3, 2))
        );
        assert_eq!(
            Number::Fraction(Rational64::new(1, 2)).mul(&Number::Real(3.0)),
            Number::Real(1.5)
        );
        assert_eq!(
            Number::Fraction(Rational64::new(1, 2))
                .mul(&Number::Complex(Complex64::new(1.0, -3.0))),
            Number::Complex(Complex64::new(0.5, -1.5))
        );
        assert_eq!(
            Number::Fraction(Rational64::new(1, 2)).mul(&Number::Fraction(Rational64::new(3, 4))),
            Number::Fraction(Rational64::new(3, 8))
        );
    }

    // Implementation for Expr
    #[test]
    fn test_impl_expr() {
        let int: Arc<dyn Expr> = Number::Integer(1).into();
        let real: Arc<dyn Expr> = Number::Real(0.0).into();
        let cmplx: Arc<dyn Expr> = Number::Complex(Complex64::new(1.0, -1.0)).into();
        let frac: Arc<dyn Expr> = Number::Fraction(Rational64::new(3, 4)).into();

        assert_eq!(int.hash_key(), "Integer(1)");
        assert_eq!(real.hash_key(), "Real(0)");
        assert_eq!(cmplx.hash_key(), "Complex(1-1i)");
        assert_eq!(frac.hash_key(), "Fraction(3/4)");

        assert!(int.is_scalar());
        assert!(real.is_scalar());
        assert!(cmplx.is_scalar());
        assert!(frac.is_scalar());

        assert!(int == Number::from_i64(1));
        assert!(real == Number::from_f64(0.0));
        assert!(cmplx == Number::from_complex(Complex64::new(1.0, -1.0)));
        assert!(frac == Number::from_rational(Rational64::new(3, 4)));

        assert!(int != real.clone());
        assert!(int != cmplx.clone());
        assert!(int != frac.clone());
        assert!(real != cmplx.clone());
        assert!(real != frac.clone());
        assert!(cmplx != frac.clone());

        assert_eq!(format!("{}", int), "1");
        assert_eq!(format!("{}", real), "0");
        assert_eq!(format!("{}", cmplx), "1 + -1i");
        assert_eq!(format!("{}", frac), "3/4");
    }

    // Test serialization and deserialization via `serde_json`
    #[test]
    fn test_serialization() {
        let n = Number::Fraction(Rational64::new(3, 7));
        let json = serde_json::to_string(&n).unwrap();
        let recovered: Number = serde_json::from_str(&json).unwrap();
        assert_eq!(n, recovered);

        let n = Number::Complex(Complex64::new(1.0, -1.0));
        let json = serde_json::to_string(&n).unwrap();
        let recovered: Number = serde_json::from_str(&json).unwrap();
        assert_eq!(n, recovered);
    }

    // Test utils: interning, downcast, type and identity check
    #[test]
    fn test_utils() {
        let int = Number::from_i64(1);
        let int1 = Number::Integer(1).into();
        let int2 = Number::from_i64(7);

        assert!(Arc::ptr_eq(&int, &int1));
        assert!(!Arc::ptr_eq(&int, &int2));

        let real = Number::from_f64(0.0);
        let real1 = Number::Real(0.0).into();
        let real2 = Number::from_f64(3.14);

        assert!(Arc::ptr_eq(&real, &real1));
        assert!(!Arc::ptr_eq(&real, &real2));

        let cmplx = Number::from_complex(Complex64::new(1.0, 0.0));
        let cmplx1 = Number::Complex(Complex64::new(1.0, 0.0)).into();
        let cmplx2 = Number::from_complex(Complex64::new(3.0, 4.0));

        assert!(Arc::ptr_eq(&cmplx, &cmplx1));
        assert!(!Arc::ptr_eq(&cmplx, &cmplx2));

        let frac = Number::from_rational(Rational64::new(3, 4));
        let frac1 = Number::Fraction(Rational64::new(3, 4)).into();
        let frac2 = Number::from_rational(Rational64::new(5, 6));

        assert!(Arc::ptr_eq(&frac, &frac1));
        assert!(!Arc::ptr_eq(&frac, &frac2));

        let mut num = downcast_from_arc::<Number>(&int).unwrap();
        assert_eq!(num, &Number::Integer(1));

        num = downcast_from_arc::<Number>(&real).unwrap();
        assert_eq!(num, &Number::Real(0.0));

        num = downcast_from_arc::<Number>(&cmplx).unwrap();
        assert_eq!(num, &Number::Complex(Complex64::new(1.0, 0.0)));

        num = downcast_from_arc::<Number>(&frac).unwrap();
        assert_eq!(num, &Number::Fraction(Rational64::new(3, 4)));

        assert!(is_expr_type::<Number>(&int));
        assert!(is_expr_type::<Number>(&real));
        assert!(is_expr_type::<Number>(&cmplx));
        assert!(is_expr_type::<Number>(&frac));

        assert!(!is_zero_expr(&int));
        assert!(is_zero_expr(&real));
        assert!(!is_zero_expr(&cmplx));
        assert!(!is_zero_expr(&frac));

        assert!(is_one_expr(&int));
        assert!(!is_one_expr(&real));
        assert!(is_one_expr(&cmplx));
        assert!(!is_one_expr(&frac));
    }
}
