use pyo3::prelude::*;

use tinned::ExpAdjointMap;

use crate::core::{errors::to_pyerr, expr::PyExpr};
use crate::perturbations::pert_multichain::PyPertMultichain;

/// Create an ExpAdjointMap expression.
///
/// Args:
///   generator: Generator expression.
///   target: Target expression.
///   generator_derivative_commute: Optional bool indicating whether the
///                                 generator and its derivatives commute.
///   left_action: Optional bool.
///                If True: exp(X)*Y*exp(-X).
///                If False: exp(-X)*Y*exp(X).
///   max_commutator_order: Optional u32 truncation.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
#[pyo3(signature = (generator, target, generator_derivative_commute=None, left_action=None, max_commutator_order=None))]
pub fn exp_adjoint_map_new(
    generator: PyExpr,
    target: PyExpr,
    generator_derivative_commute: Option<bool>,
    left_action: Option<bool>,
    max_commutator_order: Option<u32>,
) -> PyResult<PyExpr> {
    let mut b = ExpAdjointMap::builder(
        generator.inner().clone(),
        target.inner().clone(),
        generator_derivative_commute,
    );

    if let Some(v) = left_action {
        b = b.left_action(v);
    }
    if let Some(v) = max_commutator_order {
        b = b.max_commutator_order(v);
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

/// Create an ExpAdjointMap expression with target as the time evolution of generator.
///
/// Args:
///   generator: Generator expression.
///   is_forward: If True uses i*d/dt, otherwise -i*d/dt.
///   generator_derivative_commute: Optional bool indicating whether the
///                                 generator and its derivatives commute.
///   left_action: Optional bool. If True: exp(X)*Y*exp(-X), otherwise exp(-X)*Y*exp(X).
///   max_commutator_order: Optional u32 truncation.
///
/// Returns:
///   A PyExpr wrapping the constructed expression (interned).
#[pyfunction]
#[pyo3(signature = (generator, is_forward, generator_derivative_commute=None, left_action=None, max_commutator_order=None))]
pub fn exp_adjoint_map_time_evolution_new(
    generator: PyExpr,
    is_forward: bool,
    generator_derivative_commute: Option<bool>,
    left_action: Option<bool>,
    max_commutator_order: Option<u32>,
) -> PyResult<PyExpr> {
    let mut b = ExpAdjointMap::builder_time_evolution(
        generator.inner().clone(),
        is_forward,
        generator_derivative_commute,
    );

    if let Some(v) = left_action {
        b = b.left_action(v);
    }
    if let Some(v) = max_commutator_order {
        b = b.max_commutator_order(v);
    }

    let out = b.build().map_err(to_pyerr)?;
    Ok(PyExpr::new(out))
}

impl_expr_getter_interface!(
    fn_name = exp_adjoint_map_generator,
    fn_doc = impl_expr_getter_doc!("generator", ExpAdjointMap),
    expr_ty = ExpAdjointMap,
    out_ty = PyExpr,
    body = |op: &ExpAdjointMap| Ok(PyExpr::new(op.generator().clone()))
);

impl_expr_getter_interface!(
    fn_name = exp_adjoint_map_target,
    fn_doc = impl_expr_getter_doc!("target", ExpAdjointMap),
    expr_ty = ExpAdjointMap,
    out_ty = PyExpr,
    body = |op: &ExpAdjointMap| Ok(PyExpr::new(op.target().clone()))
);

impl_expr_getter_interface!(
    fn_name = exp_adjoint_map_is_time_evolution,
    fn_doc =
        impl_expr_getter_doc!("whether target is the time-differentiated generator", ExpAdjointMap),
    expr_ty = ExpAdjointMap,
    out_ty = bool,
    body = |op: &ExpAdjointMap| Ok(op.is_time_evolution())
);

impl_expr_getter_interface!(
    fn_name = exp_adjoint_map_left_action,
    fn_doc = impl_expr_getter_doc!("Boolean value of left action", ExpAdjointMap),
    expr_ty = ExpAdjointMap,
    out_ty = bool,
    body = |op: &ExpAdjointMap| Ok(op.left_action())
);

impl_expr_getter_interface!(
    fn_name = exp_adjoint_map_max_commutator_order,
    fn_doc = impl_expr_getter_doc!("maximum fold", ExpAdjointMap),
    expr_ty = ExpAdjointMap,
    out_ty = u32,
    body = |op: &ExpAdjointMap| Ok(op.max_commutator_order())
);

impl_expr_getter_interface!(
    fn_name = exp_adjoint_map_at_zero_perturbations,
    fn_doc = impl_expr_getter_doc!("whether evaluated at zero-field strength", ExpAdjointMap),
    expr_ty = ExpAdjointMap,
    out_ty = bool,
    body = |op: &ExpAdjointMap| Ok(op.at_zero_perturbations())
);

impl_expr_getter_interface!(
    fn_name = exp_adjoint_map_result,
    fn_doc = impl_expr_getter_doc!("result expression", ExpAdjointMap),
    expr_ty = ExpAdjointMap,
    out_ty = PyExpr,
    body = |op: &ExpAdjointMap| Ok(PyExpr::new(op.result().clone()))
);

impl_expr_getter_interface!(
    fn_name = exp_adjoint_map_derivative,
    fn_doc = impl_expr_getter_doc!("derivative", ExpAdjointMap),
    expr_ty = ExpAdjointMap,
    out_ty = PyPertMultichain,
    body = |op: &ExpAdjointMap| Ok(PyPertMultichain::new(op.derivative().clone()))
);

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(exp_adjoint_map_new, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_time_evolution_new, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_generator, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_target, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_is_time_evolution, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_left_action, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_max_commutator_order, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_at_zero_perturbations, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_result, m)?)?;
    m.add_function(wrap_pyfunction!(exp_adjoint_map_derivative, m)?)?;

    Ok(())
}
