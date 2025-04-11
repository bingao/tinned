use num_complex::Complex64;
use num_rational::Rational64;
use serde_json;
use std::sync::Arc;

use tinned::utils::{downcast_from_arc, is_expr_type, is_one_expr, is_zero_expr};
use tinned::{Expr, Number};

#[test]
fn test_is_zero() {
    assert!(Number::Integer(0).is_zero());
    assert!(Number::Real(0.0).is_zero());
    assert!(Number::Fraction(Rational64::new(0, 1)).is_zero());
    assert!(Number::Complex(Complex64::new(0.0, 0.0)).is_zero());
    assert!(!Number::Integer(1).is_zero());
}

#[test]
fn test_is_one() {
    assert!(Number::Integer(1).is_one());
    assert!(Number::Real(1.0).is_one());
    assert!(Number::Fraction(Rational64::new(1, 1)).is_one());
    assert!(Number::Complex(Complex64::new(1.0, 0.0)).is_one());
    assert!(!Number::Integer(0).is_one());
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
    assert_eq!(Number::Real(2.5).add(&Number::Fraction(Rational64::new(3, 4))), Number::Real(3.25));

    assert_eq!(
        Number::Complex(Complex64::new(2.0, 3.0)).add(&Number::Integer(3)),
        Number::Complex(Complex64::new(5.0, 3.0))
    );
    assert_eq!(
        Number::Complex(Complex64::new(2.0, 3.0)).add(&Number::Real(3.0)),
        Number::Complex(Complex64::new(5.0, 3.0))
    );
    assert_eq!(
        Number::Complex(Complex64::new(2.0, 3.0)).add(&Number::Complex(Complex64::new(1.0, -3.0))),
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
    assert_eq!(Number::Fraction(Rational64::new(1, 2)).add(&Number::Real(3.0)), Number::Real(3.5));
    assert_eq!(
        Number::Fraction(Rational64::new(1, 2)).add(&Number::Complex(Complex64::new(1.0, -3.0))),
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
        Number::Complex(Complex64::new(2.0, 3.0)).mul(&Number::Complex(Complex64::new(1.0, -3.0))),
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
    assert_eq!(Number::Fraction(Rational64::new(1, 2)).mul(&Number::Real(3.0)), Number::Real(1.5));
    assert_eq!(
        Number::Fraction(Rational64::new(1, 2)).mul(&Number::Complex(Complex64::new(1.0, -3.0))),
        Number::Complex(Complex64::new(0.5, -1.5))
    );
    assert_eq!(
        Number::Fraction(Rational64::new(1, 2)).mul(&Number::Fraction(Rational64::new(3, 4))),
        Number::Fraction(Rational64::new(3, 8))
    );
}

#[test]
fn test_equality() {
    assert_eq!(Number::Integer(5), Number::Integer(5));
    assert_ne!(Number::Integer(5), Number::Real(5.0)); // they are different variants
    assert_eq!(Number::Real(5.0), Number::Real(5.0));
    assert_eq!(
        Number::Complex(Complex64::new(2.0, 3.0)),
        Number::Complex(Complex64::new(2.0, 3.0))
    );
    assert_eq!(Number::Fraction(Rational64::new(1, 2)), Number::Fraction(Rational64::new(1, 2)));
}

#[test]
fn test_display() {
    assert_eq!(format!("{}", Number::Integer(42)), "42");
    assert_eq!(format!("{}", Number::Real(3.14)), "3.14");
    assert_eq!(format!("{}", Number::Fraction(Rational64::new(3, 2))), "3/2");
    assert_eq!(format!("{}", Number::Complex(Complex64::new(1.0, 2.0))), "1 + 2i");
}

#[test]
fn test_serialization_roundtrip() {
    let n = Number::Fraction(Rational64::new(3, 7));
    let json = serde_json::to_string(&n).unwrap();
    let recovered: Number = serde_json::from_str(&json).unwrap();
    assert_eq!(n, recovered);

    let n = Number::Complex(Complex64::new(1.0, -1.0));
    let json = serde_json::to_string(&n).unwrap();
    let recovered: Number = serde_json::from_str(&json).unwrap();
    assert_eq!(n, recovered);
}

#[test]
fn test_expr_trait() {
    let int: Arc<dyn Expr> = Number::Integer(1).into();
    let real: Arc<dyn Expr> = Number::Real(0.0).into();
    let cmplx: Arc<dyn Expr> = Number::Complex(Complex64::new(1.0, -1.0)).into();
    let frac: Arc<dyn Expr> = Number::Fraction(Rational64::new(3, 4)).into();

    assert!(downcast_from_arc::<Number>(&int).is_some());
    assert!(downcast_from_arc::<Number>(&real).is_some());
    assert!(downcast_from_arc::<Number>(&cmplx).is_some());
    assert!(downcast_from_arc::<Number>(&frac).is_some());

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

    assert!(!(int == real.clone()));
    assert!(!(int == cmplx.clone()));
    assert!(!(int == frac.clone()));
    assert!(!(real == cmplx.clone()));
    assert!(!(real == frac.clone()));
    assert!(!(cmplx == frac.clone()));

    assert_eq!(format!("{}", int), "1");
    assert_eq!(format!("{}", real), "0");
    assert_eq!(format!("{}", cmplx), "1 + -1i");
    assert_eq!(format!("{}", frac), "3/4");
}

#[test]
fn test_utils() {
    let int: Arc<dyn Expr> = Number::Integer(1).into();
    let real: Arc<dyn Expr> = Number::Real(0.0).into();
    let cmplx: Arc<dyn Expr> = Number::Complex(Complex64::new(1.0, 0.0)).into();
    let frac: Arc<dyn Expr> = Number::Fraction(Rational64::new(3, 4)).into();

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
