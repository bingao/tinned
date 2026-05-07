use std::sync::Arc;

use crate::core::{Expr, TinnedError};
use crate::public::{generic_expression_error, is_zero_expr};

// Applies a fallible operation to one child expression and handles the common
// rebuild logic used by unary expression types.
//
// Semantics:
// - The operation is applied to the child argument.
// - If the resulting child is zero, the owner expression collapses to the
//   zero expression supplied by `zero_expr`.
// - If the child is unchanged, the original owner expression is returned.
// - Otherwise, the owner expression is rebuilt from the transformed child via
//   `build_expr`.
#[inline]
pub(crate) fn transform_unary_any_zero(
    owner: &dyn Expr,
    arg: &Arc<dyn Expr>,
    operation: impl FnOnce(&Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
    err_message: impl Into<String>,
    build_expr: impl FnOnce(Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
    zero_expr: impl FnOnce() -> Result<Arc<dyn Expr>, TinnedError>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    let new_arg = operation(arg)
        .map_err(|e| generic_expression_error(err_message, owner, Some(Box::new(e))))?;

    if is_zero_expr(&new_arg, None) {
        zero_expr()
    } else if &new_arg == arg {
        Ok(owner.clone_expr())
    } else {
        build_expr(new_arg)
    }
}

// Applies fallible operations to two child expressions and handles the common
// rebuild logic used by binary expression types where any zero child collapses
// the whole owner expression to zero.
//
// Semantics:
// - The operation is applied to both child arguments.
// - If the resulting first child is zero, the owner expression collapses to the
//   zero expression supplied by `zero_expr` and the second child is not visited.
// - If the resulting second child is zero, the owner expression collapses to the
//   zero expression supplied by `zero_expr`.
// - If both children are unchanged, the original owner expression is returned.
// - Otherwise, the owner expression is rebuilt from the transformed children via
//   `build_expr`.
#[inline]
pub(crate) fn transform_binary_any_zero(
    owner: &dyn Expr,
    first_arg: &Arc<dyn Expr>,
    second_arg: &Arc<dyn Expr>,
    operation: impl Fn(&Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
    err_message: impl Into<String>,
    build_expr: impl FnOnce(Arc<dyn Expr>, Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
    zero_expr: impl FnOnce() -> Result<Arc<dyn Expr>, TinnedError>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    let err_message = err_message.into();

    let new_first = operation(first_arg).map_err(|e| {
        generic_expression_error(
            format!("{} for the first argument", err_message),
            owner,
            Some(Box::new(e)),
        )
    })?;

    if is_zero_expr(&new_first, None) {
        return zero_expr();
    }

    let new_second = operation(second_arg).map_err(|e| {
        generic_expression_error(
            format!("{} for the second argument", err_message),
            owner,
            Some(Box::new(e)),
        )
    })?;

    if is_zero_expr(&new_second, None) {
        return zero_expr();
    }

    if &new_first == first_arg && &new_second == second_arg {
        Ok(owner.clone_expr())
    } else {
        build_expr(new_first, new_second)
    }
}

// Applies fallible operations to two child expressions and handles the common
// rebuild logic used by binary expression types where the owner expression
// collapses to zero only when both transformed children are zero.
//
// Semantics:
// - The operation is applied to both child arguments.
// - If both transformed children are zero, the owner expression collapses to the
//   zero expression supplied by `zero_expr`.
// - If only the transformed first child is zero, the original first child is
//   retained.
// - If only the transformed second child is zero, the original second child is
//   retained.
// - If neither effective child changed, the original owner expression is
//   returned.
// - Otherwise, the owner expression is rebuilt from the effective children via
//   `build_expr`.
#[inline]
pub(crate) fn transform_binary_all_zero(
    owner: &dyn Expr,
    first_arg: &Arc<dyn Expr>,
    second_arg: &Arc<dyn Expr>,
    operation: impl Fn(&Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
    err_message: impl Into<String>,
    build_expr: impl FnOnce(Arc<dyn Expr>, Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
    zero_expr: impl FnOnce() -> Result<Arc<dyn Expr>, TinnedError>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    let err_message = err_message.into();

    let new_first = operation(first_arg).map_err(|e| {
        generic_expression_error(
            format!("{} for the first argument", err_message),
            owner,
            Some(Box::new(e)),
        )
    })?;

    let new_second = operation(second_arg).map_err(|e| {
        generic_expression_error(
            format!("{} for the second argument", err_message),
            owner,
            Some(Box::new(e)),
        )
    })?;

    let first_is_zero = is_zero_expr(&new_first, None);
    let second_is_zero = is_zero_expr(&new_second, None);

    if first_is_zero && second_is_zero {
        return zero_expr();
    }

    let effective_first = if first_is_zero {
        first_arg.clone()
    } else {
        new_first
    };

    let effective_second = if second_is_zero {
        second_arg.clone()
    } else {
        new_second
    };

    let first_changed = &effective_first != first_arg;

    let second_changed = &effective_second != second_arg;

    if !first_changed && !second_changed {
        Ok(owner.clone_expr())
    } else {
        build_expr(effective_first, effective_second)
    }
}

// Applies a fallible operation to terms of a (matrix) addition and handles the
// common rebuild logic where the owner expression collapses to zero only when
// all transformed terms are zero.
//
// Semantics:
// - The operation is applied to terms one by one.
// - If all transformed terms are zero, the owner expression collapses to the
//   zero expression supplied by `zero_expr`.
// - If no term changed, the original owner expression is returned.
// - Otherwise, the owner expression is rebuilt from the transformed terms via
//   `build_expr`.
// - Errors from `operation` are wrapped with `err_message`, argument context,
//   and the owner expression as context.
#[inline]
pub(crate) fn transform_add_all_zero(
    owner: &dyn Expr,
    terms: &[Arc<dyn Expr>],
    operation: impl Fn(&Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
    err_message: impl Into<String>,
    build_expr: impl FnOnce(Vec<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError>,
    zero_expr: impl FnOnce() -> Result<Arc<dyn Expr>, TinnedError>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    let err_message = err_message.into();

    let mut new_terms = Vec::with_capacity(terms.len());
    let mut all_zero = true;
    let mut changed = false;

    for (index, term) in terms.iter().enumerate() {
        let new_term = operation(term).map_err(|e| {
            generic_expression_error(
                format!("{err_message} for term {index}"),
                owner,
                Some(Box::new(e)),
            )
        })?;

        if is_zero_expr(&new_term, None) {
            changed = true;
        } else {
            all_zero = false;

            if &new_term != term {
                changed = true;
            }

            new_terms.push(new_term);
        }
    }

    if all_zero {
        return zero_expr();
    }

    if changed {
        build_expr(new_terms)
    } else {
        Ok(owner.clone_expr())
    }
}

// Applies a fallible operation to coefficient and factors of a (matrix)
// multiplication and handles the common rebuild logic where any zero child
// collapses the whole owner expression to zero.
//
// Semantics:
// - The operation is applied to `coefficient` first.
// - The operation is then applied to `factors` from left to right.
// - If any transformed child is zero, the owner expression collapses to the zero
//   expression supplied by `zero_expr`, and no remaining children are visited.
// - If every transformed child is unchanged, the original owner expression is
//   returned.
// - Otherwise, the owner expression is rebuilt from the transformed
//   coefficient and factors via `build_expr`.
// - Errors from `operation` are wrapped with `err_message`, argument context,
//   and the owner expression as context.
#[inline]
pub(crate) fn transform_mul_any_zero(
    owner: &dyn Expr,
    coefficient: &Arc<dyn Expr>,
    factors: &[Arc<dyn Expr>],
    operation: impl Fn(&Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
    err_message: impl Into<String>,
    build_expr: impl FnOnce(Vec<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError>,
    zero_expr: impl FnOnce() -> Result<Arc<dyn Expr>, TinnedError>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    let err_message = err_message.into();

    let new_coef = operation(coefficient).map_err(|e| {
        generic_expression_error(
            format!("{err_message} for the coefficient"),
            owner,
            Some(Box::new(e)),
        )
    })?;

    if is_zero_expr(&new_coef, None) {
        return zero_expr();
    }

    let mut new_factors = Vec::with_capacity(factors.len() + 1);
    let mut changed = &new_coef != coefficient;

    new_factors.push(new_coef);

    for (index, factor) in factors.iter().enumerate() {
        let new_factor = operation(factor).map_err(|e| {
            generic_expression_error(
                format!("{err_message} for factor {index}"),
                owner,
                Some(Box::new(e)),
            )
        })?;

        if is_zero_expr(&new_factor, None) {
            return zero_expr();
        }

        if &new_factor != factor {
            changed = true;
        }

        new_factors.push(new_factor);
    }

    if changed {
        build_expr(new_factors)
    } else {
        Ok(owner.clone_expr())
    }
}

// Applies a fallible operation to coefficient and factors of a (matrix)
// multiplication and handles the common rebuild logic where the owner
// expression collapses to zero only when all transformed children are zero.
//
// Semantics:
// - The operation is applied to `coefficient` first.
// - The operation is then applied to `factors` from left to right.
// - If all transformed children are zero, the owner expression collapses to the
//   zero expression supplied by `zero_expr`.
// - If only some transformed children are zero, those children are treated as
//   unchanged and their original child expressions are retained.
// - If no effective child changed, the original owner expression is returned.
// - Otherwise, the owner expression is rebuilt from the effective coefficient
//   and factors via `build_expr`.
// - Errors from `operation` are wrapped with `err_message`, argument context,
//   and the owner expression as context.
#[inline]
pub(crate) fn transform_mul_all_zero(
    owner: &dyn Expr,
    coefficient: &Arc<dyn Expr>,
    factors: &[Arc<dyn Expr>],
    operation: impl Fn(&Arc<dyn Expr>) -> Result<Arc<dyn Expr>, TinnedError>,
    err_message: impl Into<String>,
    build_expr: impl FnOnce(Vec<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError>,
    zero_expr: impl FnOnce() -> Result<Arc<dyn Expr>, TinnedError>,
) -> Result<Arc<dyn Expr>, TinnedError> {
    let err_message = err_message.into();

    let new_coef = operation(coefficient).map_err(|e| {
        generic_expression_error(
            format!("{err_message} for the coefficient"),
            owner,
            Some(Box::new(e)),
        )
    })?;

    let coef_is_zero = is_zero_expr(&new_coef, None);

    let effective_coef = if coef_is_zero {
        coefficient.clone()
    } else {
        new_coef
    };

    let mut effective_factors = Vec::with_capacity(factors.len());
    let mut all_zero = coef_is_zero;
    let mut changed = &effective_coef != coefficient;

    effective_factors.push(effective_coef);

    for (index, factor) in factors.iter().enumerate() {
        let new_factor = operation(factor).map_err(|e| {
            generic_expression_error(
                format!("{err_message} for factor {index}"),
                owner,
                Some(Box::new(e)),
            )
        })?;

        if is_zero_expr(&new_factor, None) {
            effective_factors.push(factor.clone());
        } else {
            all_zero = false;

            if &new_factor != factor {
                changed = true;
            }

            effective_factors.push(new_factor);
        }
    }

    if all_zero {
        return zero_expr();
    }

    if changed {
        build_expr(effective_factors)
    } else {
        Ok(owner.clone_expr())
    }
}
