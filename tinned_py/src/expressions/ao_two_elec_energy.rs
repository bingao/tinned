use pyo3::prelude::*;

use tinned::{AoTwoElecEnergy, AoTwoElecMatrix, expression_error};

use crate::core::{errors::to_pyerr, expr::PyExpr};
use crate::perturbations::pert_multichain::PyPertMultichain;

/// Create a AoTwoElecEnergy expression.
///
/// Args:
///   name: Energy name.
///   inner_density: Inner density expression (must be WfnParameter or ResidueParameter; ZeroOperator yields Number(0)).
///   outer_density: Optional outer density expression. Defaults to inner_density.
///   allow_density_swap: Optional bool controlling equality behavior. Defaults to True.
///   dependencies: Optional dependencies. Defaults to empty.
///   derivative: Optional derivative. Defaults to empty.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned), or Number(0) when invalid/zero per rules.
#[pyfunction]
#[pyo3(signature = (
    name,
    inner_density,
    outer_density=None,
    allow_density_swap=None,
    dependencies=None,
    derivative=None
))]
pub fn ao_two_elec_energy_new(
    name: String,
    inner_density: PyExpr,
    outer_density: Option<PyExpr>,
    allow_density_swap: Option<bool>,
    dependencies: Option<&PyPertMultichain>,
    derivative: Option<&PyPertMultichain>,
) -> PyResult<PyExpr> {
    let mut b = AoTwoElecEnergy::builder(name, inner_density.inner().clone());

    if let Some(v) = outer_density {
        b = b.outer_density(v.inner().clone());
    }
    if let Some(v) = allow_density_swap {
        b = b.allow_density_swap(v);
    }
    if let Some(v) = dependencies {
        b = b.dependencies(v.inner().clone());
    }
    if let Some(v) = derivative {
        b = b.derivative(v.inner().clone());
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Create an AoTwoElecEnergy expression from an AoTwoElecMatrix expression.
///
/// Args:
///   two_elec_op: A AoTwoElecMatrix expression.
///   outer_density: Optional outer density expression. Defaults to operator density.
///   allow_density_swap: Optional bool controlling equality behavior. Defaults to True.
///   dependencies: Optional dependencies. Defaults to operator dependencies.
///   derivative: Optional derivative. Defaults to operator derivative.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned), or Number(0) when invalid/zero per rules.
#[pyfunction]
#[pyo3(signature = (
    two_elec_op,
    outer_density=None,
    allow_density_swap=None,
    dependencies=None,
    derivative=None
))]
pub fn ao_two_elec_energy_from_operator(
    two_elec_op: PyExpr,
    outer_density: Option<PyExpr>,
    allow_density_swap: Option<bool>,
    dependencies: Option<&PyPertMultichain>,
    derivative: Option<&PyPertMultichain>,
) -> PyResult<PyExpr> {
    let inner = two_elec_op.inner().clone();

    let op_ref = inner.as_any().downcast_ref::<AoTwoElecMatrix>().ok_or_else(|| {
        to_pyerr(expression_error(
            "ao_two_elec_energy_from_operator() expected a AoTwoElecMatrix expression",
            &inner,
            None,
        ))
    })?;

    let mut b = AoTwoElecEnergy::builder_from_operator(op_ref);

    if let Some(v) = outer_density {
        b = b.outer_density(v.inner().clone());
    }
    if let Some(v) = allow_density_swap {
        b = b.allow_density_swap(v);
    }
    if let Some(v) = dependencies {
        b = b.dependencies(v.inner().clone());
    }
    if let Some(v) = derivative {
        b = b.derivative(v.inner().clone());
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = ao_two_elec_energy_name,
    fn_doc = impl_expr_getter_doc!("name", AoTwoElecEnergy),
    expr_ty = AoTwoElecEnergy,
    out_ty = String,
    body = |op: &AoTwoElecEnergy| Ok(op.name().to_string())
);

impl_expr_getter_interface!(
    fn_name = ao_two_elec_energy_inner_density,
    fn_doc = impl_expr_getter_doc!("inner density", AoTwoElecEnergy),
    expr_ty = AoTwoElecEnergy,
    out_ty = PyExpr,
    body = |op: &AoTwoElecEnergy| Ok(PyExpr::new(op.inner_density().clone()))
);

impl_expr_getter_interface!(
    fn_name = ao_two_elec_energy_outer_density,
    fn_doc = impl_expr_getter_doc!("outer density", AoTwoElecEnergy),
    expr_ty = AoTwoElecEnergy,
    out_ty = PyExpr,
    body = |op: &AoTwoElecEnergy| Ok(PyExpr::new(op.outer_density().clone()))
);

impl_expr_getter_interface!(
    fn_name = ao_two_elec_energy_allow_density_swap,
    fn_doc = impl_expr_getter_doc!("whether density swapping is allowed", AoTwoElecEnergy),
    expr_ty = AoTwoElecEnergy,
    out_ty = bool,
    body = |op: &AoTwoElecEnergy| Ok(op.allow_density_swap())
);

impl_expr_getter_interface!(
    fn_name = ao_two_elec_energy_dependencies,
    fn_doc = impl_expr_getter_doc!("dependencies", AoTwoElecEnergy),
    expr_ty = AoTwoElecEnergy,
    out_ty = PyPertMultichain,
    body = |op: &AoTwoElecEnergy| Ok(PyPertMultichain::new(op.dependencies().clone()))
);

impl_expr_getter_interface!(
    fn_name = ao_two_elec_energy_derivative,
    fn_doc = impl_expr_getter_doc!("derivative", AoTwoElecEnergy),
    expr_ty = AoTwoElecEnergy,
    out_ty = PyPertMultichain,
    body = |op: &AoTwoElecEnergy| Ok(PyPertMultichain::new(op.derivative().clone()))
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ao_two_elec_energy_new, m)?)?;
    m.add_function(wrap_pyfunction!(ao_two_elec_energy_from_operator, m)?)?;
    m.add_function(wrap_pyfunction!(ao_two_elec_energy_name, m)?)?;
    m.add_function(wrap_pyfunction!(ao_two_elec_energy_inner_density, m)?)?;
    m.add_function(wrap_pyfunction!(ao_two_elec_energy_outer_density, m)?)?;
    m.add_function(wrap_pyfunction!(ao_two_elec_energy_allow_density_swap, m)?)?;
    m.add_function(wrap_pyfunction!(ao_two_elec_energy_dependencies, m)?)?;
    m.add_function(wrap_pyfunction!(ao_two_elec_energy_derivative, m)?)?;

    Ok(())
}
