use pyo3::basic::CompareOp;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use tinned::Perturbation;

use crate::core::expr::PyExpr;

#[pyclass(module = "tinned", name = "Perturbation", from_py_object)]
#[derive(Clone)]
pub struct PyPerturbation {
    inner: Arc<Perturbation>,
}

impl PyPerturbation {
    pub fn new(inner: Arc<Perturbation>) -> Self {
        Self {
            inner,
        }
    }

    pub fn inner(&self) -> &Arc<Perturbation> {
        &self.inner
    }
}

#[pymethods]
impl PyPerturbation {
    #[new]
    fn __new__(name: String, frequency: PyExpr) -> Self {
        Self {
            inner: Perturbation::new(name, frequency.inner().clone()),
        }
    }

    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    #[getter]
    fn frequency(&self) -> PyExpr {
        PyExpr::new(self.inner.frequency().clone())
    }

    fn __repr__(&self) -> String {
        format!("Perturbation({})", self.inner)
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __hash__(&self) -> isize {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        self.inner().as_ref().hash(&mut h);
        h.finish() as isize
    }

    fn __richcmp__(&self, other: PyRef<'_, PyPerturbation>, op: CompareOp) -> PyResult<bool> {
        match op {
            CompareOp::Eq => Ok(self.inner.as_ref() == other.inner.as_ref()),
            CompareOp::Ne => Ok(self.inner.as_ref() != other.inner.as_ref()),
            _ => Err(PyTypeError::new_err("Perturbation only supports == and !=")),
        }
    }

    fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(&*self.inner).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[staticmethod]
    fn from_json(s: &str) -> PyResult<Self> {
        let p: Perturbation =
            serde_json::from_str(s).map_err(|e| PyValueError::new_err(e.to_string()))?;

        // Re-intern by going through Perturbation::new (recommended),
        // but we need frequency as Arc<dyn Expr>. If your deserialization already
        // returns that correctly, you can rewrap directly. Otherwise, adapt.
        //
        // If deserializing Perturbation yields name + frequency in the right types:
        let arc = Perturbation::new(p.name().to_string(), p.frequency().clone());
        Ok(Self::new(arc))
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPerturbation>()?;
    Ok(())
}
