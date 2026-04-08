use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use pyo3::Py;
use pyo3::basic::CompareOp;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use tinned::{Expr, Perturbation, expr_from_json, expr_to_json};

use crate::core::errors::to_pyerr;
use crate::perturbations::perturbation::PyPerturbation;
use crate::public::number_tolerance::PyNumberTolerance;

#[pyclass(module = "tinned", name = "Expr", from_py_object)]
#[derive(Clone)]
pub struct PyExpr {
    inner: Arc<dyn Expr>,
}

impl PyExpr {
    pub fn new(inner: Arc<dyn Expr>) -> Self {
        Self {
            inner,
        }
    }

    pub fn inner(&self) -> &Arc<dyn Expr> {
        &self.inner
    }
}

#[pymethods]
impl PyExpr {
    fn __repr__(&self) -> String {
        format!("Expr({})", self.inner)
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn is_scalar(&self) -> bool {
        self.inner.is_scalar()
    }

    fn has_unperturbed_term(&self) -> bool {
        self.inner.has_unperturbed_term()
    }

    /// Performs a conditional canonicalization to zero, such as setting
    /// `TimeEvolution` and unperturbed `BasisTimeEvolution` to zero, and
    /// undifferentiated perturbing operators to zero.
    ///
    /// Args:
    ///   freq_tol: Optional NumberTolerance.
    ///
    /// Returns:
    ///   A PyExpr wrapping the canonicalized expression.
    #[pyo3(signature = (freq_tol=None))]
    fn substitute_zero_perturbations(&self, freq_tol: Option<PyNumberTolerance>) -> PyResult<PyExpr> {
        let tol = freq_tol.map(|t| t.into_inner());
        let out = self.inner.substitute_zero_perturbations(tol).map_err(to_pyerr)?;
        Ok(PyExpr::new(out))
    }

    /// Differentiate with respect to a Perturbation.
    ///
    /// Args:
    ///   s: Perturbation.
    ///
    /// Returns:
    ///   A PyExpr wrapping the differentiated expression.
    fn differentiate(&self, s: &PyPerturbation) -> PyResult<PyExpr> {
        let out = self.inner.differentiate(s.inner()).map_err(to_pyerr)?;
        Ok(PyExpr::new(out))
    }

    /// Eliminate derivatives of a response parameter from the expression.
    ///
    /// Args:
    ///   parameter: A LagMultiplier or WfnParameter expression.
    ///   perturbations: Iterable of Perturbation.
    ///   min_order: Minimum order.
    ///
    /// Returns:
    ///   A PyExpr wrapping the eliminated expression.
    fn eliminate(
        &self,
        parameter: PyExpr,
        perturbations: Vec<PyPerturbation>,
        min_order: u32,
    ) -> PyResult<PyExpr> {
        let perts: Vec<Arc<Perturbation>> =
            perturbations.into_iter().map(|p| p.inner().clone()).collect();

        let out = self.inner.eliminate(parameter.inner(), &perts, min_order).map_err(to_pyerr)?;

        Ok(PyExpr::new(out))
    }

    /// Check if any expression in set exists in the current expression.
    ///
    /// Args:
    ///   set: Iterable of Expr.
    ///   include_derivatives: Whether considering derivatives (including order
    ///                        0) of expressions in set when checking existence.
    ///
    /// Returns:
    ///   True if any exists, otherwise False.
    fn exist_any(&self, set: Vec<PyExpr>, include_derivatives: bool) -> PyResult<bool> {
        let rust_set: HashSet<Arc<dyn Expr>> = set.into_iter().map(|e| e.inner().clone()).collect();
        Ok(self.inner.exist_any(&rust_set, include_derivatives))
    }

    /// Find the given expression and all its higher-order superchain matches.
    ///
    /// Args:
    ///   s: Expr to match.
    ///
    /// Returns:
    ///   A dict mapping total_order (int) to a list of matching Expr.
    fn find_superchains<'py>(&self, py: Python<'py>, s: PyExpr) -> PyResult<Py<PyDict>> {
        let m: BTreeMap<u32, HashSet<Arc<dyn Expr>>> = self.inner.find_superchains(s.inner());

        let out = PyDict::new(py);
        for (order, set) in m {
            let list: Vec<PyExpr> = set.into_iter().map(PyExpr::new).collect();
            out.set_item(order, list)?;
        }

        Ok(out.into())
    }

    /// Remove all expressions in set from the current expression.
    ///
    /// Args:
    ///   set: Iterable of Expr.
    ///
    /// Returns:
    ///   A PyExpr wrapping the resulting expression.
    fn remove(&self, set: Vec<PyExpr>) -> PyResult<PyExpr> {
        let rust_set: HashSet<Arc<dyn Expr>> = set.into_iter().map(|e| e.inner().clone()).collect();
        let out = self.inner.remove(&rust_set).map_err(to_pyerr)?;
        Ok(PyExpr::new(out))
    }

    /// Replace expressions according to a mapping.
    ///
    /// Args:
    ///   map: dict[Expr, Expr]
    ///   include_derivatives: Whether considering derivatives (including order
    ///                        0) of expressions in set when checking existence.
    ///
    /// Returns:
    ///   A PyExpr wrapping the replaced expression.
    fn replace(&self, map: &Bound<'_, PyDict>, include_derivatives: bool) -> PyResult<PyExpr> {
        let mut rust_map: HashMap<Arc<dyn Expr>, Arc<dyn Expr>> = HashMap::new();

        for (k, v) in map.iter() {
            let key: PyExpr = k.extract()?;
            let val: PyExpr = v.extract()?;
            rust_map.insert(key.inner().clone(), val.inner().clone());
        }

        let out = self.inner.replace(&rust_map, include_derivatives).map_err(to_pyerr)?;
        Ok(PyExpr::new(out))
    }

    /// Retain expressions in set while remove others from the current expression.
    ///
    /// Args:
    ///   set: Iterable of Expr.
    ///   include_derivatives: Whether considering derivatives (including order
    ///                        0) of expressions in set when checking existence.
    ///
    /// Returns:
    ///   A PyExpr wrapping the retained expression.
    fn retain(&self, set: Vec<PyExpr>, include_derivatives: bool) -> PyResult<PyExpr> {
        let rust_set: HashSet<Arc<dyn Expr>> = set.into_iter().map(|e| e.inner().clone()).collect();
        let out = self.inner.retain(&rust_set, include_derivatives).map_err(to_pyerr)?;
        Ok(PyExpr::new(out))
    }

    fn __hash__(&self) -> isize {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        self.inner().hash(&mut h); // uses dyn Expr Hash impl (hash_key)
        h.finish() as isize
    }

    fn __richcmp__(&self, other: PyRef<'_, PyExpr>, op: CompareOp) -> PyResult<bool> {
        let a = self.inner();
        let b = other.inner();

        match op {
            CompareOp::Eq => Ok(a.as_ref() == b.as_ref()),
            CompareOp::Ne => Ok(a.as_ref() != b.as_ref()),
            CompareOp::Lt => Ok(a.as_ref().cmp(b.as_ref()) == Ordering::Less),
            CompareOp::Le => Ok(a.as_ref().cmp(b.as_ref()) != Ordering::Greater),
            CompareOp::Gt => Ok(a.as_ref().cmp(b.as_ref()) == Ordering::Greater),
            CompareOp::Ge => Ok(a.as_ref().cmp(b.as_ref()) != Ordering::Less),
        }
    }

    fn to_json(&self) -> PyResult<String> {
        expr_to_json(self.inner()).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[staticmethod]
    fn from_json(s: &str) -> PyResult<Self> {
        let expr = expr_from_json(s).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(PyExpr::new(expr))
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyExpr>()?;
    Ok(())
}
