use num_complex::Complex64;
use num_rational::Rational64;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyAny;

use tinned::Number;

use crate::core::expr::PyExpr;

fn rational64_from_py(obj: &Bound<'_, PyAny>) -> PyResult<Rational64> {
    // Accept:
    // 1) fractions.Fraction
    // 2) Any object with numerator and denominator attributes
    // 3) A 2-tuple (numerator, denominator)

    if let Ok((n, d)) = obj.extract::<(i64, i64)>() {
        if d == 0 {
            return Err(PyValueError::new_err("Denominator must be non-zero"));
        }
        return Ok(Rational64::new(n, d));
    }

    let has_num = obj.hasattr("numerator")?;
    let has_den = obj.hasattr("denominator")?;
    if has_num && has_den {
        let n = obj.getattr("numerator")?.extract::<i64>()?;
        let d = obj.getattr("denominator")?.extract::<i64>()?;
        if d == 0 {
            return Err(PyValueError::new_err("Denominator must be non-zero"));
        }
        return Ok(Rational64::new(n, d));
    }

    Err(PyTypeError::new_err(
        "Expected a rational as fractions.Fraction, (numerator, denominator), or an object with numerator/denominator",
    ))
}

/// Create a symbolic integer expression.
///
/// Agrs:
///   n: Integer value.
///
/// Returns
///   A PyExpr wrapping the constructed Number (interned).
#[pyfunction]
pub fn number_from_i64(n: i64) -> PyResult<PyExpr> {
    Ok(PyExpr::new(Number::from_i64(n)))
}

/// Create a symbolic real number expression.
///
/// Agrs:
///   x: Real value.
///
/// Returns
///   A PyExpr wrapping the constructed Number (interned).
#[pyfunction]
pub fn number_from_f64(x: f64) -> PyResult<PyExpr> {
    Ok(PyExpr::new(Number::from_f64(x)))
}

/// Create a symbolic complex number expression.
///
/// Agrs:
///   z: Complex number.
///
/// Returns
///   A PyExpr wrapping the constructed Number (interned).
#[pyfunction]
pub fn number_from_complex(z: Complex64) -> PyResult<PyExpr> {
    // PyO3 can extract Python complex into num_complex::Complex64 if the feature is enabled.
    Ok(PyExpr::new(Number::from_complex(z)))
}

/// Create a symbolic rational number expression.
///
/// Accepted inputs:
/// - fractions.Fraction
/// - (numerator, denominator) tuple
/// - Any object with attributes 'numerator' and 'denominator'
///
/// Agrs:
///   r: Rational value.
///
/// Returns
///   A PyExpr wrapping the constructed Number (interned).
#[pyfunction]
pub fn number_from_rational(r: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    let rat = rational64_from_py(r)?;
    Ok(PyExpr::new(Number::from_rational(rat)))
}

/// Return the symbolic zero.
///
/// Returns
///   A PyExpr wrapping the constructed Number (interned).
#[pyfunction]
pub fn number_zero() -> PyResult<PyExpr> {
    Ok(PyExpr::new(Number::zero()))
}

/// Return the symbolic one.
///
/// Returns
///   A PyExpr wrapping the constructed Number (interned).
#[pyfunction]
pub fn number_one() -> PyResult<PyExpr> {
    Ok(PyExpr::new(Number::one()))
}

/// Return the symbolic minus one.
///
/// Returns
///   A PyExpr wrapping the constructed Number (interned).
#[pyfunction]
pub fn number_minus_one() -> PyResult<PyExpr> {
    Ok(PyExpr::new(Number::minus_one()))
}

/// Return the symbolic one half (1/2).
///
/// Returns
///   A PyExpr wrapping the constructed Number (interned).
#[pyfunction]
pub fn number_one_half() -> PyResult<PyExpr> {
    Ok(PyExpr::new(Number::one_half()))
}

/// Return the symbolic imaginary unit i.
///
/// Returns
///   A PyExpr wrapping the constructed Number (interned).
#[pyfunction]
pub fn number_imaginary_unit() -> PyResult<PyExpr> {
    Ok(PyExpr::new(Number::imaginary_unit()))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(number_from_i64, m)?)?;
    m.add_function(wrap_pyfunction!(number_from_f64, m)?)?;
    m.add_function(wrap_pyfunction!(number_from_complex, m)?)?;
    m.add_function(wrap_pyfunction!(number_from_rational, m)?)?;

    m.add_function(wrap_pyfunction!(number_zero, m)?)?;
    m.add_function(wrap_pyfunction!(number_one, m)?)?;
    m.add_function(wrap_pyfunction!(number_minus_one, m)?)?;
    m.add_function(wrap_pyfunction!(number_one_half, m)?)?;
    m.add_function(wrap_pyfunction!(number_imaginary_unit, m)?)?;

    Ok(())
}
