macro_rules! expr_arc_ty {
    () => {
        ::std::sync::Arc<dyn $crate::core::Expr>
    };
}

macro_rules! expr_arc_ref_ty {
    () => {
        &expr_arc_ty!()
    };
}

macro_rules! expr_result_ty {
    () => {
        ::std::result::Result<expr_arc_ty!(), $crate::core::TinnedError>
    };
}

macro_rules! expr_map_ty {
    () => {
        ::std::collections::HashMap<expr_arc_ty!(), expr_arc_ty!()>
    };
}

macro_rules! expr_set_ty {
    () => {
        ::std::collections::HashSet<expr_arc_ty!()>
    };
}

macro_rules! expr_differentiation_map_ty {
    () => {
        ::std::collections::BTreeMap<u32, expr_set_ty!()>
    };
}
