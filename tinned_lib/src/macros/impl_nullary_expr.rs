macro_rules! impl_nullary_expr_type {
    ($type_name:ident, $builder_name:ident, $has_deps:tt, $is_scalar:tt) => {
        impl_nullary_expr_type!(@nullary_def_oper $type_name, $has_deps);

        impl $type_name {
            impl_nullary_expr_type!(@nullary_oper_methods $builder_name, $has_deps);

            #[inline]
            pub fn name(&self) -> &str {
                &self.name
            }

            #[inline]
            pub fn derivative(&self) -> &$crate::perturbations::PertMultichain {
                &self.derivative
            }
        }

        impl_nullary_expr_type!(@nullary_def_builder $builder_name, $has_deps);

        impl $builder_name {
            #[inline]
            pub fn derivative(
                mut self,
                derivative: $crate::perturbations::PertMultichain,
            ) -> Self {
                self.derivative = derivative;
                self
            }

            impl_nullary_expr_type!(@nullary_builder_methods $type_name, $has_deps, $is_scalar);
        }
    };

    (@nullary_def_oper $type_name:ident, true) => {
        #[derive(Clone, Debug, PartialEq, Eq, ::serde::Serialize, ::serde::Deserialize)]
        pub struct $type_name {
            name: ::std::string::String,
            dependencies: $crate::perturbations::PertMultichain,
            derivative: $crate::perturbations::PertMultichain,
            is_perturbing: bool,
        }
    };

    (@nullary_def_oper $type_name:ident, false) => {
        #[derive(Clone, Debug, PartialEq, Eq, ::serde::Serialize, ::serde::Deserialize)]
        pub struct $type_name {
            name: ::std::string::String,
            derivative: $crate::perturbations::PertMultichain,
        }
    };

    (@nullary_def_builder $builder_name:ident, true) => {
        #[derive(Debug)]
        pub struct $builder_name {
            name: ::std::string::String,
            dependencies: $crate::perturbations::PertMultichain,
            derivative: $crate::perturbations::PertMultichain,
            is_perturbing: bool,
        }
    };

    (@nullary_def_builder $builder_name:ident, false) => {
        #[derive(Debug)]
        pub struct $builder_name {
            name: ::std::string::String,
            derivative: $crate::perturbations::PertMultichain,
        }
    };

    (@nullary_oper_methods $builder_name:ident, true) => {
        #[inline]
        pub fn builder(
            name: impl ::std::convert::Into<::std::string::String>,
        ) -> $builder_name {
            $builder_name {
                name: name.into(),
                dependencies: $crate::perturbations::PertMultichain::new(),
                derivative: $crate::perturbations::PertMultichain::new(),
                is_perturbing: false,
            }
        }

        #[inline]
        fn with_derivative(
            &self,
            derivative: $crate::perturbations::PertMultichain,
        ) -> $builder_name {
            $builder_name {
                name: self.name.clone(),
                dependencies: self.dependencies.clone(),
                derivative,
                is_perturbing: self.is_perturbing,
            }
        }

        #[inline]
        pub fn dependencies(&self) -> &$crate::perturbations::PertMultichain {
            &self.dependencies
        }

        #[inline]
        pub fn is_perturbing(&self) -> bool {
            self.is_perturbing
        }
    };

    (@nullary_oper_methods $builder_name:ident, false) => {
        #[inline]
        pub fn builder(
            name: impl ::std::convert::Into<::std::string::String>,
        ) -> $builder_name {
            $builder_name {
                name: name.into(),
                derivative: $crate::perturbations::PertMultichain::new(),
            }
        }

        #[inline]
        fn with_derivative(
            &self,
            derivative: $crate::perturbations::PertMultichain,
        ) -> $builder_name {
            $builder_name {
                name: self.name.clone(),
                derivative,
            }
        }
    };

    (@nullary_builder_methods $type_name:ident, true, true) => {
        #[inline]
        pub fn dependencies(
            mut self,
            deps: $crate::perturbations::PertMultichain,
        ) -> Self {
            self.dependencies = deps;
            self
        }

        #[inline]
        pub fn is_perturbing(mut self, is_perturbing: bool) -> Self {
            self.is_perturbing = is_perturbing;
            self
        }

        #[inline]
        pub fn build(self) -> expr_result_ty!() {
            if self.dependencies.is_subchain(&self.derivative) {
                Ok($crate::internal::intern_expr(::std::sync::Arc::new($type_name {
                    name: self.name,
                    dependencies: self.dependencies,
                    derivative: self.derivative,
                    is_perturbing: self.is_perturbing,
                })))
            } else {
                Ok($crate::expressions::Number::zero())
            }
        }
    };

    (@nullary_builder_methods $type_name:ident, true, false) => {
        #[inline]
        pub fn dependencies(
            mut self,
            deps: $crate::perturbations::PertMultichain,
        ) -> Self {
            self.dependencies = deps;
            self
        }

        #[inline]
        pub fn is_perturbing(mut self, is_perturbing: bool) -> Self {
            self.is_perturbing = is_perturbing;
            self
        }

        #[inline]
        pub fn build(self) -> expr_result_ty!() {
            if self.dependencies.is_subchain(&self.derivative) {
                Ok($crate::internal::intern_expr(::std::sync::Arc::new($type_name {
                    name: self.name,
                    dependencies: self.dependencies,
                    derivative: self.derivative,
                    is_perturbing: self.is_perturbing,
                })))
            } else {
                Ok($crate::expressions::ZeroOperator::new())
            }
        }
    };

    (@nullary_builder_methods $type_name:ident, false, true) => {
        compile_error!(
            "impl_nullary_expr_type!(...) does not support has_deps = false and is_scalar = true"
        );
    };

    (@nullary_builder_methods $type_name:ident, false, false) => {
        #[inline]
        pub fn build(self) -> expr_result_ty!() {
            Ok($crate::internal::intern_expr(::std::sync::Arc::new($type_name {
                name: self.name,
                derivative: self.derivative,
            })))
        }
    };
}

macro_rules! impl_nullary_expr_traits {
    ($type_name:ident, $has_deps:tt, $is_scalar:tt) => {
        impl $crate::core::ExprInternal for $type_name {
            impl_expr_internal_methods!($type_name, true);

            impl_nullary_expr_traits!(@nullary_hash_key $type_name, $has_deps);

            #[inline]
            fn total_order(&self) -> u32 {
                self.derivative.total_order()
            }

            #[inline]
            fn deep_eq_superchains(
                &self,
                other: &expr_arc_ty!(),
            ) -> bool {
                if let Some(op) = $crate::public::downcast_from_arc::<$type_name>(other) {
                    impl_nullary_expr_traits!(
                        @nullary_deep_eq_superchains self,
                        op,
                        $has_deps
                    )
                } else {
                    false
                }
            }

            #[inline]
            fn eq_by_superchains(
                &self,
                other: &expr_arc_ty!(),
            ) -> bool {
                self.deep_eq_superchains(other)
            }

            #[inline]
            fn replace_expr_fields(
                &self,
                _map: &expr_map_ty!(),
                _exact_equality: bool,
            ) -> expr_result_ty!() {
                Ok(self.clone_expr())
            }

            #[inline]
            fn retain_expr_fields(
                &self,
                _expr: &expr_arc_ty!(),
                _exact_equality: bool,
            ) -> expr_result_ty!() {
                impl_zero_expr!($is_scalar)
            }
        }

        #[::typetag::serde]
        impl $crate::core::Expr for $type_name {
            impl_nullary_expr_common_methods!($type_name, $is_scalar);

            impl_nullary_expr_traits!(@nullary_apply_zero_rules $type_name, $has_deps, $is_scalar);

            #[inline]
            fn differentiate(
                &self,
                s: &::std::sync::Arc<$crate::perturbations::Perturbation>,
            ) -> expr_result_ty!() {
                let new_deriv = self.derivative.with_added_perturbation(s);
                self.with_derivative(new_deriv).build()
            }

            impl_nullary_expr_traits!(@nullary_eliminate $type_name, $has_deps);
        }

        impl_nullary_expr_traits!(@nullary_display $type_name, $has_deps);
    };

    (@nullary_hash_key $type_name:ident, true) => {
        #[inline]
        fn hash_key(&self) -> ::std::string::String {
            ::std::format!(
                "{}({}; [{}]; [{}]; {})",
                stringify!($type_name),
                self.name,
                self.dependencies.hash_key(),
                self.derivative.hash_key(),
                self.is_perturbing,
            )
        }
    };

    (@nullary_hash_key $type_name:ident, false) => {
        #[inline]
        fn hash_key(&self) -> ::std::string::String {
            ::std::format!(
                "{}({}; [{}])",
                stringify!($type_name),
                self.name,
                self.derivative.hash_key(),
            )
        }
    };

    (@nullary_deep_eq_superchains $self:ident, $op:ident, true) => {
        $self.name == $op.name
            && $self.dependencies == $op.dependencies
            && $self.derivative.is_subchain(&$op.derivative)
            && $self.is_perturbing == $op.is_perturbing
    };

    (@nullary_deep_eq_superchains $self:ident, $op:ident, false) => {
        $self.name == $op.name
            && $self.derivative.is_subchain(&$op.derivative)
    };

    (@nullary_apply_zero_rules $type_name:ident, true, $is_scalar:tt) => {
        #[inline]
        fn apply_zero_rules(
            &self,
            _freq_tol: ::std::option::Option<$crate::public::NumberTolerance>,
        ) -> expr_result_ty!() {
            // For a perturbing operator, it is non-zero only when it is
            // differentiated with at least one dependency
            if self.is_perturbing && self.derivative.is_empty() {
                impl_zero_expr!($is_scalar)
            } else {
                Ok(self.clone_expr())
            }
        }
    };

    (@nullary_apply_zero_rules $type_name:ident, false, true) => {};

    (@nullary_apply_zero_rules $type_name:ident, false, false) => {};

    (@nullary_eliminate $type_name:ident, true) => {};

    (@nullary_eliminate $type_name:ident, false) => {
        #[inline]
        fn eliminate(
            &self,
            parameter: &expr_arc_ty!(),
            perturbations: &[::std::sync::Arc<$crate::perturbations::Perturbation>],
            min_order: u32,
        ) -> expr_result_ty!() {
            if let Some(op) = $crate::public::downcast_from_arc::<$type_name>(parameter) {
                if self.name == op.name {
                    let map = self.derivative.get_map_clone();

                    let order: u32 = perturbations
                        .iter()
                        .map(|p| *map.get(p).unwrap_or(&0))
                        .sum();

                    if order >= min_order
                        && order <= perturbations.len() as u32
                    {
                        return Ok($crate::expressions::ZeroOperator::new());
                    }
                }
            }

            Ok(self.clone_expr())
        }
    };

    (@nullary_display $type_name:ident, true) => {
        impl ::std::fmt::Display for $type_name {
            fn fmt(
                &self,
                f: &mut ::std::fmt::Formatter,
            ) -> ::std::fmt::Result {
                if self.derivative.is_empty() {
                    ::std::write!(f, "{}({})", self.name, self.is_perturbing)
                } else {
                    ::std::write!(
                        f,
                        "{}({})^({})",
                        self.name,
                        self.is_perturbing,
                        self.derivative
                    )
                }
            }
        }
    };

    (@nullary_display $type_name:ident, false) => {
        impl ::std::fmt::Display for $type_name {
            fn fmt(
                &self,
                f: &mut ::std::fmt::Formatter,
            ) -> ::std::fmt::Result {
                if self.derivative.is_empty() {
                    ::std::write!(f, "{}", self.name)
                } else {
                    ::std::write!(
                        f,
                        "{}^({})",
                        self.name,
                        self.derivative
                    )
                }
            }
        }
    };
}

macro_rules! impl_nullary_expr_common_methods {
    ($type_name:ident, $is_scalar:tt) => {
        impl_expr_common_methods!($is_scalar);

        #[inline]
        fn remove(&self, set: &expr_set_ty!()) -> expr_result_ty!() {
            if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                impl_zero_expr!($is_scalar)
            } else {
                Ok(self.clone_expr())
            }
        }
    };
}
