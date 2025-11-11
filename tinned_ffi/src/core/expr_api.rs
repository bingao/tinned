use safer_ffi::prelude::*;
use std::collections::{BTreeMap, HashSet};
use std::sync::Arc;

use tinned::core::{Expr, TinnedError};
use tinned::public::{NumberTolerance, generic_error};

use crate::c_support::{
    tinned_string_from_cstr, tinned_string_to_cstr, try_from_handle, try_with_handle,
};
use crate::core::{
    ExprBox, ExprHandle, ExprSlice, TinnedErrorBox, expr_map_from_slices, expr_set_from_slice,
    tinned_error_new,
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
fn ffi_return_string(
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
fn ffi_return_exprbox(
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
    ffi_return_string(h, "tinned_expr_type_name", out_err, |expr| Ok(expr.type_name().to_string()))
}

#[ffi_export]
pub fn tinned_expr_hash_key(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    ffi_return_string(h, "tinned_expr_hash_key", out_err, |expr| Ok(expr.hash_key()))
}

#[ffi_export]
pub fn tinned_expr_display(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    ffi_return_string(h, "tinned_expr_display", out_err, |expr| Ok(format!("{}", expr)))
}

#[ffi_export]
pub fn tinned_expr_serialize_json(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    ffi_return_string(h, "tinned_expr_serialize_json", out_err, |expr| {
        serde_json::to_string(expr.as_ref()).map_err(|err| {
            generic_error("Failed to serialize expression to JSON", Some(Box::new(err)))
        })
    })
}

// Whether the expression is scalar. Returns false on error. NULL input.
#[ffi_export]
pub fn tinned_expr_is_scalar(
    h: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    match try_with_handle(h, "tinned_expr_is_scalar", "ExprHandle", |eh| {
        let expr = eh.as_ref();
        Ok(expr.is_scalar())
    }) {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            false
        },
    }
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

// Cleans `TemporumOperator` and unperturbed `TemporumOverlap` objects.
#[ffi_export]
pub fn tinned_expr_clean_temporum(
    h: Option<&ExprHandle>,
    tol: Option<&NumberToleranceHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let tol_opt: Option<NumberTolerance> = tol.map(|t| t.as_ref().clone());
    ffi_return_exprbox(h, "tinned_expr_clean_temporum", out_err, move |expr| {
        expr.clean_temporum(tol_opt)
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
    ffi_return_exprbox(h, "tinned_expr_differentiate", out_err, move |expr| {
        expr.differentiate(&pert)
    })
}

// Eliminates a given response `parameter`'s derivatives from the expression.
#[ffi_export]
pub fn tinned_expr_eliminate(
    h: Option<&ExprHandle>,
    parameter: Option<&ExprHandle>,
    perturbations: Option<PerturbationSlice<'_>>,
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

    match try_with_handle(h, "tinned_expr_eliminate", "ExprHandle", |eh| {
        eh.as_ref().eliminate(&param, &perts, min_order).map(|e| ExprBox::new(ExprHandle::new(e)))
    }) {
        Ok(b) => Some(b),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Checks if any expression in `set` exists in the current expression.
#[ffi_export]
pub fn tinned_expr_exist_any(
    h: Option<&ExprHandle>,
    set: Option<ExprSlice<'_>>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> bool {
    let set_hs = match set {
        Some(slice) => match expr_set_from_slice(slice, "tinned_expr_exist_any") {
            Ok(s) => s,
            Err(e) => {
                tinned_error_new(out_err, e);
                return false;
            },
        },
        None => HashSet::new(),
    };

    match try_with_handle(h, "tinned_expr_exist_any", "ExprHandle", |eh| {
        Ok(eh.as_ref().exist_any(&set_hs))
    }) {
        Ok(v) => v,
        Err(e) => {
            tinned_error_new(out_err, e);
            false
        },
    }
}

// Finds a given expression `s` and all its higher-order "differentiated" ones in the current expression.
#[ffi_export]
pub fn tinned_expr_find_superchains_json(
    h: Option<&ExprHandle>,
    s: Option<&ExprHandle>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<char_p::Box> {
    let s_expr = match try_from_handle(s, "tinned_expr_find_superchains_json", "ExprHandle", |eh| {
        eh.clone_arc()
    }) {
        Ok(x) => x,
        Err(e) => {
            tinned_error_new(out_err, e);
            return None;
        },
    };

    match try_with_handle(h, "tinned_expr_find_superchains_json", "ExprHandle", |eh| {
        let map = eh.as_ref().find_superchains(&s_expr);

        // Convert HashSet<Arc<dyn Expr>> -> Vec<Arc<dyn Expr>> for JSON.
        let mut out: BTreeMap<u32, Vec<Arc<dyn Expr>>> = BTreeMap::new();
        for (k, vset) in map {
            let mut v: Vec<Arc<dyn Expr>> = vset.into_iter().collect();
            // Optional: stable order by hash_key to keep output deterministic.
            v.sort_by(|a, b| a.hash_key().cmp(&b.hash_key()));
            out.insert(k, v);
        }

        serde_json::to_string(&out).map(tinned_string_to_cstr).map_err(|err| {
            generic_error("Failed to serialize find_superchains map to JSON", Some(Box::new(err)))
        })
    }) {
        Ok(s) => Some(s),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

// Removes all expressions in `set` from the current expression.
#[ffi_export]
pub fn tinned_expr_remove(
    h: Option<&ExprHandle>,
    set: Option<ExprSlice<'_>>,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let set_hs = match set {
        Some(slc) => match expr_set_from_slice(slc, "tinned_expr_remove") {
            Ok(s) => s,
            Err(e) => {
                tinned_error_new(out_err, e);
                return None;
            },
        },
        None => Default::default(),
    };
    ffi_return_exprbox(h, "tinned_expr_remove", out_err, move |expr| expr.remove(&set_hs))
}

#[ffi_export]
pub fn tinned_expr_replace(
    h: Option<&ExprHandle>,
    keys: Option<ExprSlice<'_>>,
    values: Option<ExprSlice<'_>>,
    exact_equality: bool,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let map = match (keys, values) {
        (Some(k), Some(v)) => match expr_map_from_slices(k, v, "tinned_expr_replace") {
            Ok(m) => m,
            Err(e) => {
                tinned_error_new(out_err, e);
                return None;
            },
        },
        _ => {
            tinned_error_new(
                out_err,
                generic_error("tinned_expr_replace: keys/values must be non-NULL", None),
            );
            return None;
        },
    };

    match try_with_handle(h, "tinned_expr_replace", "ExprHandle", |eh| {
        eh.as_ref().replace(&map, exact_equality).map(|e| ExprBox::new(ExprHandle::new(e)))
    }) {
        Ok(b) => Some(b),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
}

#[ffi_export]
pub fn tinned_expr_retain(
    h: Option<&ExprHandle>,
    set: Option<ExprSlice<'_>>,
    exact_equality: bool,
    out_err: Option<Out<'_, TinnedErrorBox>>,
) -> Option<ExprBox> {
    let set_hs = match set {
        Some(slice) => match expr_set_from_slice(slice, "tinned_expr_retain") {
            Ok(s) => s,
            Err(e) => {
                tinned_error_new(out_err, e);
                return None;
            },
        },
        None => HashSet::new(),
    };

    match try_with_handle(h, "tinned_expr_retain", "ExprHandle", |eh| {
        eh.as_ref().retain(&set_hs, exact_equality).map(|e| ExprBox::new(ExprHandle::new(e)))
    }) {
        Ok(b) => Some(b),
        Err(e) => {
            tinned_error_new(out_err, e);
            None
        },
    }
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
