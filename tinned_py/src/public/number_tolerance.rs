use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

use tinned::number_tolerance::{
    NumberTolerance, get_number_tolerance as rs_get_number_tolerance,
    set_number_tolerance as rs_set_number_tolerance,
};

#[pyclass(name = "NumberTolerance", frozen, from_py_object)]
#[derive(Clone)]
pub struct PyNumberTolerance {
    inner: NumberTolerance,
}

impl PyNumberTolerance {
    #[inline]
    pub fn inner(&self) -> &NumberTolerance {
        &self.inner
    }

    #[inline]
    pub fn into_inner(self) -> NumberTolerance {
        self.inner
    }
}

#[pymethods]
impl PyNumberTolerance {
    #[new]
    #[pyo3(signature = (abs_error, rel_error))]
    pub fn new(abs_error: f64, rel_error: f64) -> PyResult<Self> {
        if abs_error < 0.0 {
            return Err(PyValueError::new_err("NumberTolerance abs_error must be non-negative"));
        }
        if rel_error < 0.0 {
            return Err(PyValueError::new_err("NumberTolerance rel_error must be non-negative"));
        }

        Ok(Self {
            inner: NumberTolerance::new(abs_error, rel_error),
        })
    }

    #[staticmethod]
    pub fn zero() -> Self {
        Self {
            inner: NumberTolerance::zero(),
        }
    }

    #[getter]
    pub fn abs_error(&self) -> f64 {
        self.inner.abs_error()
    }

    #[getter]
    pub fn rel_error(&self) -> f64 {
        self.inner.rel_error()
    }

    pub fn max_abs_error(&self, a: f64, b: f64) -> f64 {
        self.inner.max_abs_error(a, b)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "NumberTolerance(abs_error={}, rel_error={})",
            self.inner.abs_error(),
            self.inner.rel_error()
        )
    }

    pub fn __richcmp__(&self, other: PyRef<PyNumberTolerance>, op: pyo3::basic::CompareOp) -> bool {
        match op {
            pyo3::basic::CompareOp::Eq => {
                self.inner.abs_error().to_bits() == other.inner.abs_error().to_bits()
                    && self.inner.rel_error().to_bits() == other.inner.rel_error().to_bits()
            },
            pyo3::basic::CompareOp::Ne => {
                self.inner.abs_error().to_bits() != other.inner.abs_error().to_bits()
                    || self.inner.rel_error().to_bits() != other.inner.rel_error().to_bits()
            },
            _ => false,
        }
    }
}

#[pyfunction(name = "get_number_tolerance")]
pub fn get_number_tolerance_py() -> PyResult<PyNumberTolerance> {
    // rs_get_number_tolerance() panics only if the lock is poisoned.
    // Convert that to a Python exception.
    let tol = std::panic::catch_unwind(|| rs_get_number_tolerance())
        .map_err(|_| PyRuntimeError::new_err("Failed to read global NUMBER_TOLERANCE"))?;

    Ok(PyNumberTolerance {
        inner: tol,
    })
}

#[pyfunction(name = "set_number_tolerance")]
pub fn set_number_tolerance_py(new_tol: PyNumberTolerance) -> PyResult<()> {
    let tol = new_tol.inner().clone();

    let result = std::panic::catch_unwind(|| rs_set_number_tolerance(tol))
        .map_err(|_| PyRuntimeError::new_err("Failed to write global NUMBER_TOLERANCE"))?;

    let _ = result;
    Ok(())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyNumberTolerance>()?;
    m.add_function(wrap_pyfunction!(get_number_tolerance_py, m)?)?;
    m.add_function(wrap_pyfunction!(set_number_tolerance_py, m)?)?;

    Ok(())
}
