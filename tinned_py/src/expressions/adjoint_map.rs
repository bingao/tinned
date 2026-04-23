use pyo3::conversion::{FromPyObject, IntoPyObject};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyString};

use tinned::{AdjointMap, AdjointMode};

use crate::core::{errors::to_pyerr, expr::PyExpr};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PyAdjointMode(AdjointMode);

impl From<PyAdjointMode> for AdjointMode {
    fn from(value: PyAdjointMode) -> Self {
        value.0
    }
}

impl From<AdjointMode> for PyAdjointMode {
    fn from(value: AdjointMode) -> Self {
        Self(value)
    }
}

impl<'py> FromPyObject<'_, 'py> for PyAdjointMode {
    type Error = PyErr;

    fn extract(obj: Borrowed<'_, 'py, PyAny>) -> Result<Self, Self::Error> {
        let s = obj.extract::<&str>()?;

        let mode = match s {
            "commutative" => AdjointMode::Commutative,
            "symmetrized" => AdjointMode::Symmetrized,
            "ordered" => AdjointMode::Ordered,
            _ => {
                return Err(PyValueError::new_err(
                    "adjoint_mode must be 'commutative', 'symmetrized', or 'ordered'",
                ));
            },
        };

        Ok(Self(mode))
    }
}

impl<'py> IntoPyObject<'py> for PyAdjointMode {
    type Target = PyString;
    type Output = Bound<'py, PyString>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        let s = match self.0 {
            AdjointMode::Commutative => "commutative",
            AdjointMode::Symmetrized => "symmetrized",
            AdjointMode::Ordered => "ordered",
        };

        Ok(PyString::new(py, s))
    }
}

/// Create an AdjointMap expression.
///
/// Args:
///   generators: A list of non-scalar operator expressions.
///   target: A non-scalar operator expression.
///   left_action: If provided, sets whether the adjoint map acts
///                from the left (default True).
///   adjoint_mode: mode of an adjoint map or its generators which can be
///                 either "commutative", "symmetrized", or "ordered" (default
///                 "commutative").
///
/// Returns:
///   A PyExpr wrapping the constructed AdjointMap (interned).
#[pyfunction]
#[pyo3(signature = (generators, target, left_action=true, adjoint_mode=PyAdjointMode(AdjointMode::Commutative)))]
pub fn adjoint_map_new(
    generators: Vec<PyExpr>,
    target: &PyExpr,
    left_action: bool,
    adjoint_mode: PyAdjointMode,
) -> PyResult<PyExpr> {
    let rust_generators = generators.into_iter().map(|t| t.inner().clone()).collect();
    let rust_target = target.inner().clone();
    let adjoint_mode: AdjointMode = adjoint_mode.into();

    let out = AdjointMap::new(rust_generators, rust_target, Some(left_action), Some(adjoint_mode))
        .map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = adjoint_map_generators,
    fn_doc = impl_expr_getter_doc!("generators", AdjointMap),
    expr_ty = AdjointMap,
    out_ty = Vec<PyExpr>,
    body = |op: &AdjointMap| Ok(op.generators().iter().cloned().map(PyExpr::new).collect())
);

impl_expr_getter_interface!(
    fn_name = adjoint_map_adjoint_mode,
    fn_doc = impl_expr_getter_doc!("Mode of the adjoint map", AdjointMap),
    expr_ty = AdjointMap,
    out_ty = PyAdjointMode,
    body = |op: &AdjointMap| Ok(PyAdjointMode(op.adjoint_mode()))
);

impl_expr_getter_interface!(
    fn_name = adjoint_map_target,
    fn_doc = impl_expr_getter_doc!("target", AdjointMap),
    expr_ty = AdjointMap,
    out_ty = PyExpr,
    body = |op: &AdjointMap| Ok(PyExpr::new(op.target().clone()))
);

impl_expr_getter_interface!(
    fn_name = adjoint_map_left_action,
    fn_doc = impl_expr_getter_doc!("Boolean value of left action", AdjointMap),
    expr_ty = AdjointMap,
    out_ty = bool,
    body = |op: &AdjointMap| Ok(op.left_action())
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(adjoint_map_new, m)?)?;
    m.add_function(wrap_pyfunction!(adjoint_map_generators, m)?)?;
    m.add_function(wrap_pyfunction!(adjoint_map_adjoint_mode, m)?)?;
    m.add_function(wrap_pyfunction!(adjoint_map_target, m)?)?;
    m.add_function(wrap_pyfunction!(adjoint_map_left_action, m)?)?;
    Ok(())
}
