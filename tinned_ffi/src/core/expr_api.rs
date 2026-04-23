use safer_ffi::prelude::*;
use std::collections::HashSet;
use std::sync::Arc;

use tinned::core::{Expr, TinnedError};
use tinned::public::{NumberTolerance, generic_error};

use crate::c_support::{
    tinned_string_from_cstr, tinned_string_to_cstr, try_from_handle, try_with_handle,
};
use crate::core::{
    ExprBox, ExprHandle, ExprSlice, ExprSuperchainBox, ExprSuperchainHandle, TinnedErrorBox,
    expr_map_from_slices, expr_set_from_slice, tinned_error_new,
};
use crate::perturbations::{PerturbationHandle, PerturbationSlice, perturbation_vec_from_slice};
use crate::public::NumberToleranceHandle;

#[inline]
fn with_expr_arc<R>(
    h: Option<&ExprHandle>,
    caller: &'static str,
    f: impl FnOnce(Arc<dyn Expr>) -> Result<R, TinnedError>,
) -> Result<R, TinnedError> {
    try_with_handle(h, caller, "ExprHandle", |eh| {
        let expr = eh.clone_arc();
        f(expr)
    })
}

#[inline]
fn ffi_expr_return_val<R>(
    h: Option<&ExprHandle>,
    caller: &'static str,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    f: impl FnOnce(Arc<dyn Expr>) -> Result<R, TinnedError>,
) -> R
where
    R: Default,
{
    match with_expr_arc(h, caller, f) {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            R::default()
        },
    }
}

#[inline]
fn ffi_expr_return_cstr(
    h: Option<&ExprHandle>,
    caller: &'static str,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    f: impl FnOnce(Arc<dyn Expr>) -> Result<String, TinnedError>,
) -> Option<char_p::Box> {
    match with_expr_arc(h, caller, f).map(tinned_string_to_cstr) {
        Ok(s) => Some(s),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

#[inline]
fn ffi_expr_return_exprbox(
    h: Option<&ExprHandle>,
    caller: &'static str,
    out_err: Option<Out<'_, TinnedErrorBox>>,
    f: impl FnOnce(Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
) -> Option<ExprBox> {
    match with_expr_arc(h, caller, f).map(|arc| ExprBox::new(ExprHandle::new(arc))) {
        Ok(b) => Some(b),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

#[ffi_export]
pub fn tinned_expr_type_name(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    ffi_expr_return_cstr(h, "tinned_expr_type_name", out_err, |expr| {
        let full = expr.type_name();
        let s = full.strip_prefix("dyn ").unwrap_or(full);
        let no_generics = s.split('<').next().unwrap_or(s);
        Ok(no_generics.rsplit("::").next().unwrap_or(no_generics).to_string())
    })
}

#[ffi_export]
pub fn tinned_expr_hash_key(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    ffi_expr_return_cstr(h, "tinned_expr_hash_key", out_err, |expr| Ok(expr.hash_key()))
}

#[ffi_export]
pub extern "C" fn tinned_expr_eq(
    lhs: Option<&ExprHandle>,
    rhs: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    ffi_expr_return_val(lhs, "tinned_expr_eq(lhs)", out_err, |l| {
        with_expr_arc(rhs, "tinned_expr_eq(rhs)", |r| Ok(l.as_ref() == r.as_ref()))
    })
}

#[ffi_export]
pub fn tinned_expr_display(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    ffi_expr_return_cstr(h, "tinned_expr_display", out_err, |expr| Ok(format!("{}", expr)))
}

#[ffi_export]
pub fn tinned_expr_serialize_json(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    ffi_expr_return_cstr(h, "tinned_expr_serialize_json", out_err, |expr| {
        serde_json::to_string(expr.as_ref()).map_err(|err| {
            generic_error("Failed to serialize expression to JSON", Some(Box::new(err)))
        })
    })
}

// Clone an expression (like Arc clone). Returns NULL on error / NULL input.
#[ffi_export]
pub fn tinned_expr_clone(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    match try_from_handle(h, "tinned_expr_clone", "ExprHandle", |eh| {
        ExprBox::new(ExprHandle::new(eh.clone_arc()))
    }) {
        Ok(b) => Some(b),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Whether the expression is scalar. Returns false on error. NULL input.
#[ffi_export]
pub fn tinned_expr_is_scalar(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    ffi_expr_return_val(h, "tinned_expr_is_scalar", out_err, |e| Ok(e.is_scalar()))
}

// Whether the expression has unperturbed term. Returns false on error. NULL input.
#[ffi_export]
pub fn tinned_expr_has_unperturbed_term(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    ffi_expr_return_val(h, "tinned_expr_has_unperturbed_term", out_err, |e| {
        Ok(e.has_unperturbed_term())
    })
}

// Cleans `TimeEvolution` and unperturbed `BasisTimeEvolution` objects.
#[ffi_export]
pub fn tinned_expr_substitute_zero_perturbations(
    h: Option<&ExprHandle>,
    tol: Option<&NumberToleranceHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let tol_opt: Option<NumberTolerance> = tol.map(|t| t.as_ref().clone());
    ffi_expr_return_exprbox(h, "tinned_expr_substitute_zero_perturbations", out_err, move |expr| {
        expr.substitute_zero_perturbations(tol_opt)
    })
}

// Differentiates with respect to a `Perturbation`.
#[ffi_export]
pub fn tinned_expr_differentiate(
    h: Option<&ExprHandle>,
    s: Option<&PerturbationHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let pert = match try_with_handle(s, "tinned_expr_differentiate", "PerturbationHandle", |ph| {
        Ok(ph.clone_arc())
    }) {
        Ok(p) => p,
        Err(e) => {
            tinned_error_new(out_err, e);
            return None;
        },
    };

    ffi_expr_return_exprbox(h, "tinned_expr_differentiate", out_err, move |expr| {
        expr.differentiate(pert)
    })
}

// Eliminates a given response `parameter`'s derivatives from the expression.
#[ffi_export]
pub fn tinned_expr_eliminate(
    h: Option<&ExprHandle>,
    parameter: Option<&ExprHandle>,
    perturbations: Option<&PerturbationSlice>,
    min_order: u32,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let param = match try_from_handle(parameter, "tinned_expr_eliminate", "ExprHandle", |eh| {
        eh.clone_arc()
    }) {
        Ok(p) => p,
        Err(e) => {
            tinned_error_new(out_err, e);
            return None;
        },
    };

    let perts = match perturbations {
        Some(slice) => match perturbation_vec_from_slice(slice, "tinned_expr_eliminate") {
            Ok(v) => v,
            Err(e) => {
                tinned_error_new(out_err, e);
                return None;
            },
        },
        None => Vec::new(),
    };

    ffi_expr_return_exprbox(h, "tinned_expr_eliminate", out_err, move |expr| {
        expr.eliminate(param, &perts, min_order)
    })
}

// Checks if any expression in `set` exists in the current expression.
#[ffi_export]
pub fn tinned_expr_match_any(
    h: Option<&ExprHandle>,
    set: Option<&ExprSlice>,
    include_derivatives: bool,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    let expr_set = match set {
        Some(slice) => {
            match expr_set_from_slice::<HashSet<Arc<dyn Expr>>>(slice, "tinned_expr_match_any") {
                Ok(s) => s,
                Err(e) => {
                    tinned_error_new(out_err, e);
                    return false;
                },
            }
        },
        None => HashSet::new(),
    };

    ffi_expr_return_val(h, "tinned_expr_match_any", out_err, |e| {
        Ok(e.match_any(&expr_set, include_derivatives))
    })
}

// Finds a given expression `s` and all its higher-order "differentiated" ones in the current expression.
#[ffi_export]
pub fn tinned_expr_find_all(
    h: Option<&ExprHandle>,
    s: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprSuperchainBox> {
    let s_expr =
        match try_from_handle(s, "tinned_expr_find_all_new", "ExprHandle", |eh| eh.clone_arc()) {
            Ok(x) => x,
            Err(e) => {
                tinned_error_new(out_err, e);
                return None;
            },
        };

    match try_with_handle(h, "tinned_expr_find_all_new", "ExprHandle", |eh| {
        let superchains = eh.as_ref().find_all(&s_expr);

        let mut orders: Vec<u32> = superchains.keys().copied().collect();
        orders.sort_unstable();

        let mut order_exprs = Vec::with_capacity(orders.len());
        for order in &orders {
            let mut v: Vec<Arc<dyn Expr>> =
                superchains.get(order).unwrap().iter().cloned().collect();
            v.sort_by(|a, b| a.hash_key().cmp(&b.hash_key()));
            order_exprs.push(v);
        }

        Ok(ExprSuperchainBox::new(ExprSuperchainHandle::new(orders, order_exprs)))
    }) {
        Ok(b) => Some(b),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Removes a given expressions `s` from the current expression.
#[ffi_export]
pub fn tinned_expr_remove_one(
    h: Option<&ExprHandle>,
    s: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let s_expr =
        match try_from_handle(s, "tinned_expr_remove_one", "ExprHandle", |eh| eh.clone_arc()) {
            Ok(p) => p,
            Err(e) => {
                tinned_error_new(out_err, e);
                return None;
            },
        };

    ffi_expr_return_exprbox(h, "tinned_expr_remove_one", out_err, move |expr| {
        expr.remove_one(&s_expr)
    })
}

// Removes all expressions in `set` from the current expression.
#[ffi_export]
pub fn tinned_expr_remove_all(
    h: Option<&ExprHandle>,
    set: Option<&ExprSlice>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let expr_set = match set {
        Some(slice) => {
            match expr_set_from_slice::<HashSet<Arc<dyn Expr>>>(slice, "tinned_expr_remove_all") {
                Ok(s) => s,
                Err(e) => {
                    tinned_error_new(out_err, e);
                    return None;
                },
            }
        },
        None => Default::default(),
    };

    ffi_expr_return_exprbox(h, "tinned_expr_remove_all", out_err, move |expr| {
        expr.remove_all(&expr_set)
    })
}

#[ffi_export]
pub fn tinned_expr_replace_all(
    h: Option<&ExprHandle>,
    keys: Option<&ExprSlice>,
    values: Option<&ExprSlice>,
    include_derivatives: bool,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let expr_map = match (keys, values) {
        (Some(k), Some(v)) => match expr_map_from_slices(k, v, "tinned_expr_replace_all") {
            Ok(m) => m,
            Err(e) => {
                tinned_error_new(out_err, e);
                return None;
            },
        },
        _ => {
            tinned_error_new(
                out_err,
                generic_error("tinned_expr_replace_all: keys/values must be non-NULL", None),
            );
            return None;
        },
    };

    ffi_expr_return_exprbox(h, "tinned_expr_replace_all", out_err, move |expr| {
        expr.replace_all(&expr_map, include_derivatives)
    })
}

#[ffi_export]
pub fn tinned_expr_retain_all(
    h: Option<&ExprHandle>,
    set: Option<&ExprSlice>,
    include_derivatives: bool,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let expr_set = match set {
        Some(slice) => {
            match expr_set_from_slice::<HashSet<Arc<dyn Expr>>>(slice, "tinned_expr_retain_all") {
                Ok(s) => s,
                Err(e) => {
                    tinned_error_new(out_err, e);
                    return None;
                },
            }
        },
        None => HashSet::new(),
    };

    ffi_expr_return_exprbox(h, "tinned_expr_retain_all", out_err, move |expr| {
        expr.retain_all(&expr_set, include_derivatives)
    })
}

// Deserialize from JSON. `json` is nullable `const char*`.
// Returns NULL on error. Use `tinned_expr_free` to free the result.
#[ffi_export]
pub fn tinned_expr_deserialize_json(
    json: Option<char_p::Ref<'_>>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let s = match tinned_string_from_cstr(json) {
        Some(s) => s,
        None => {
            tinned_error_new(
                out_err,
                generic_error("tinned_expr_deserialize_json: NULL or non-UTF-8 json", None),
            );
            return None;
        },
    };

    match serde_json::from_str::<Arc<dyn Expr>>(&s) {
        Ok(expr) => Some(ExprBox::new(ExprHandle::new(expr))),
        Err(err) => {
            tinned_error_new(
                out_err,
                generic_error("Failed to deserialize expression from JSON", Some(Box::new(err))),
            );
            None
        },
    }
}
