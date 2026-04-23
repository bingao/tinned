use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use std::hash::{Hash, Hasher};

use tinned::PertMultichain;

use crate::perturbations::perturbation::PyPerturbation;

#[pyclass(module = "tinned", name = "PertMultichain", from_py_object)]
#[derive(Clone)]
pub struct PyPertMultichain {
    inner: PertMultichain,
}

impl PyPertMultichain {
    #[inline]
    pub fn new(inner: PertMultichain) -> Self {
        Self {
            inner,
        }
    }

    pub fn inner(&self) -> &PertMultichain {
        &self.inner
    }
}

#[pymethods]
impl PyPertMultichain {
    #[new]
    fn __new__() -> Self {
        Self {
            inner: PertMultichain::new(),
        }
    }

    // Python: PertMultichain.from_iter([p1, p2, p1])
    #[staticmethod]
    fn from_iter(iterable: Vec<PyPerturbation>) -> PyResult<Self> {
        let perts = iterable.into_iter().map(|p| p.inner().clone()).collect::<Vec<_>>();

        Ok(Self {
            inner: PertMultichain::from_slice(&perts),
        })
    }

    fn insert(&mut self, p: &PyPerturbation) {
        self.inner.insert(p.inner().clone());
    }

    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn get_order(&self, p: &PyPerturbation) -> u32 {
        self.inner.get_order(p.inner())
    }

    fn total_order(&self) -> u32 {
        self.inner.total_order()
    }

    // Return owned Python objects: Vec[ Perturbation ]
    // Add `py: Python` parameter instead of Python::with_gil.
    fn keys(&self, py: Python<'_>) -> PyResult<Vec<Py<PyPerturbation>>> {
        self.inner.keys().into_iter().map(|arc_p| Py::new(py, PyPerturbation::new(arc_p))).collect()
    }

    fn to_list(&self, py: Python<'_>) -> PyResult<Vec<Py<PyPerturbation>>> {
        self.inner
            .to_vec()
            .into_iter()
            .map(|arc_p| Py::new(py, PyPerturbation::new(arc_p)))
            .collect()
    }

    fn is_subchain(&self, sub: &PyPertMultichain) -> bool {
        self.inner.is_subchain(&sub.inner)
    }

    fn is_superchain(&self, sup: &PyPertMultichain) -> bool {
        self.inner.is_superchain(&sup.inner)
    }

    fn has_overlap(&self, other: &PyPertMultichain) -> bool {
        self.inner.has_overlap(&other.inner)
    }

    fn complement(
        &self,
        py: Python<'_>,
        other: &PyPertMultichain,
    ) -> PyResult<Vec<Py<PyPerturbation>>> {
        self.inner
            .complement(&other.inner)
            .into_iter()
            .map(|arc_p| Py::new(py, PyPerturbation::new(arc_p)))
            .collect()
    }

    fn __len__(&self) -> usize {
        self.total_order() as usize
    }

    fn __repr__(&self) -> String {
        format!("PertMultichain({})", self.inner)
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __hash__(&self) -> isize {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        self.inner.hash_key().hash(&mut h);
        h.finish() as isize
    }

    fn __richcmp__(
        &self,
        other: PyRef<'_, PyPertMultichain>,
        op: pyo3::basic::CompareOp,
    ) -> PyResult<bool> {
        match op {
            pyo3::basic::CompareOp::Eq => Ok(self.inner == other.inner),
            pyo3::basic::CompareOp::Ne => Ok(self.inner != other.inner),
            _ => Err(PyTypeError::new_err("PertMultichain only supports == and !=")),
        }
    }

    fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(&self.inner()).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[staticmethod]
    fn from_json(s: &str) -> PyResult<Self> {
        let chain: PertMultichain =
            serde_json::from_str(s).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self {
            inner: chain,
        })
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPertMultichain>()?;
    Ok(())
}
