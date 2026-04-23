macro_rules! py_expr_ty {
    () => {
        $crate::core::expr::PyExpr
    };
}

macro_rules! py_expr_ref_ty {
    () => {
        &py_expr_ty!()
    };
}

macro_rules! py_pert_ty {
    () => {
        $crate::perturbations::perturbation::PyPerturbation
    };
}

macro_rules! py_pert_multichain_ty {
    () => {
        $crate::perturbations::pert_multichain::PyPertMultichain
    };
}

macro_rules! py_pert_multichain_ref_ty {
    () => {
        &py_pert_multichain_ty!()
    };
}
