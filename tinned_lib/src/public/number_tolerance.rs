use once_cell::sync::Lazy;
use std::sync::RwLock;

/// Global tolerance settings for comparing floating point numbers in `Number`.
#[derive(Clone, Debug)]
pub struct NumberTolerance {
    abs_error: f64,
    rel_error: f64,
}

impl NumberTolerance {
    #[inline]
    pub fn new(abs_error: f64, rel_error: f64) -> Self {
        if abs_error < 0.0 {
            panic!("NumberTolerance absolute error must be non-negative");
        }
        if rel_error < 0.0 {
            panic!("NumberTolerance relative error must be non-negative");
        }
        Self {
            abs_error,
            rel_error,
        }
    }

    #[inline]
    pub fn zero() -> Self {
        Self {
            abs_error: 0.0,
            rel_error: 0.0,
        }
    }

    #[inline]
    pub fn abs_error(&self) -> f64 {
        self.abs_error
    }

    #[inline]
    pub fn rel_error(&self) -> f64 {
        self.rel_error
    }

    #[inline]
    pub fn max_abs_error(&self, a: f64, b: f64) -> f64 {
        self.abs_error.max(a.abs().max(b.abs()) * self.rel_error)
    }
}

// Default values
static NUMBER_TOLERANCE: Lazy<RwLock<NumberTolerance>> =
    Lazy::new(|| RwLock::new(NumberTolerance::new(0.0, 1.0e-12)));

/// Helper to get a *copy* of current tolerance (thread-safe)
#[inline]
pub fn get_number_tolerance() -> NumberTolerance {
    NUMBER_TOLERANCE.read().expect("Failed to acquire read lock on NUMBER_TOLERANCE").clone()
}

/// Helper to update tolerance safely
#[inline]
pub fn set_number_tolerance(new_tol: NumberTolerance) {
    let mut tol =
        NUMBER_TOLERANCE.write().expect("Failed to acquire write lock on NUMBER_TOLERANCE");

    *tol = new_tol;
}
