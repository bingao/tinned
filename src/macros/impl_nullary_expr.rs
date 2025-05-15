macro_rules! impl_nullary_expr_type {
    ($type_name:ident, $builder_name:ident, $has_deps:tt, $is_scalar:tt) => {
        impl_nullary_expr_type!(@def_oper_struct $type_name, $has_deps);

        impl $type_name {
            impl_nullary_expr_type!(@impl_oper_methods $builder_name, $has_deps);

            #[inline]
            pub fn name(&self) -> &str {
                &self.name
            }

            #[inline]
            pub fn derivative(&self) -> &PertMultichain {
                &self.derivative
            }
        }

        impl_nullary_expr_type!(@def_builder_struct $builder_name, $has_deps);

        impl $builder_name {
            #[inline]
            pub fn derivative(mut self, deriv: PertMultichain) -> Self {
                self.derivative = deriv;
                self
            }

            impl_nullary_expr_type!(@impl_builder_methods $type_name, $has_deps, $is_scalar);
        }
    };

    (@def_oper_struct $type_name:ident, true) => {
        #[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
        pub struct $type_name {
            name: String,
            dependencies: PertMultichain,
            derivative: PertMultichain,
        }
    };

    (@def_oper_struct $type_name:ident, false) => {
        #[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
        pub struct $type_name {
            name: String,
            derivative: PertMultichain,
        }
    };

    (@def_builder_struct $builder_name:ident, true) => {
        #[derive(Debug)]
        pub struct $builder_name {
            name: String,
            dependencies: PertMultichain,
            derivative: PertMultichain,
        }
    };

    (@def_builder_struct $builder_name:ident, false) => {
        #[derive(Debug)]
        pub struct $builder_name {
            name: String,
            derivative: PertMultichain,
        }
    };

    (@impl_oper_methods $builder_name:ident, true) => {
        #[inline]
        pub fn builder(name: impl Into<String>) -> $builder_name {
            $builder_name {
                name: name.into(),
                dependencies: PertMultichain::new(),
                derivative: PertMultichain::new(),
            }
        }

        #[inline]
        fn builder_from(&self, derivative: PertMultichain) -> $builder_name {
            $builder_name {
                name: self.name.clone(),
                dependencies: self.dependencies.clone(),
                derivative,
            }
        }

        #[inline]
        pub fn dependencies(&self) -> &PertMultichain {
            &self.dependencies
        }
    };

    (@impl_oper_methods $builder_name:ident, false) => {
        #[inline]
        pub fn builder(name: impl Into<String>) -> $builder_name {
            $builder_name {
                name: name.into(),
                derivative: PertMultichain::new(),
            }
        }

        #[inline]
        fn builder_from(&self, derivative: PertMultichain) -> $builder_name {
            $builder_name {
                name: self.name.clone(),
                derivative,
            }
        }
    };

    (@impl_builder_methods $type_name:ident, true, true) => {
        #[inline]
        pub fn dependencies(mut self, deps: PertMultichain) -> Self {
            self.dependencies = deps;
            self
        }

        #[inline]
        pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
            if self.dependencies.is_subchain(&self.derivative) {
                Ok(crate::internal::intern_expr(Arc::new($type_name {
                    name: self.name,
                    dependencies: self.dependencies,
                    derivative: self.derivative,
                })))
            } else {
                Ok(Number::zero())
            }
        }
    };

    (@impl_builder_methods $type_name:ident, true, false) => {
        #[inline]
        pub fn dependencies(mut self, deps: PertMultichain) -> Self {
            self.dependencies = deps;
            self
        }

        #[inline]
        pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
            if self.dependencies.is_subchain(&self.derivative) {
                Ok(crate::internal::intern_expr(Arc::new($type_name {
                    name: self.name,
                    dependencies: self.dependencies,
                    derivative: self.derivative,
                })))
            } else {
                Ok(ZeroOperator::new())
            }
        }
    };

    (@impl_builder_methods $type_name:ident, false, true) => {
        compile_error!("impl_nullary_expr_type!(...) does not support has_deps = false and is_scalar = true");
    };

    (@impl_builder_methods $type_name:ident, false, false) => {
        #[inline]
        pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
            Ok(crate::internal::intern_expr(Arc::new($type_name {
                name: self.name,
                derivative: self.derivative,
            })))
        }
    };
}

macro_rules! impl_nullary_expr_traits {
    ($type_name:ident, $has_deps:tt, $is_scalar:tt) => {
        impl ExprInternal for $type_name {
            impl_expr_internal_methods!($type_name);

            impl_nullary_expr_traits!(@impl_hash_key $type_name, $has_deps);

            #[inline]
            fn total_order(&self) -> u32 {
                self.derivative.total_order()
            }

            #[inline]
            fn match_for_find_all(&self, other: &Arc<dyn Expr>) -> bool {
                if let Some(op) = downcast_from_arc::<$type_name>(other) {
                    impl_nullary_expr_traits!(@impl_match_for_find_all self, op, $has_deps)
                } else {
                    false
                }
            }

            #[inline]
            fn match_for_replace_all(&self, other: &Arc<dyn Expr>) -> bool {
                self.match_for_find_all(other)
            }
        }

        #[typetag::serde]
        impl Expr for $type_name {
            impl_nullary_expr_common_methods!($is_scalar);

            #[inline]
            fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError>
            {
                let new_deriv = self.derivative.clone_with_insert(s);

                self.builder_from(new_deriv).build()
            }

            impl_nullary_expr_traits!(@impl_eliminate $type_name, $has_deps);
        }

        impl std::fmt::Display for $type_name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                if self.derivative.is_empty() {
                    write!(f, "{}", self.name)
                } else {
                    write!(f, "{}^({})", self.name, self.derivative)
                }
            }
        }
    };

    (@impl_hash_key $type_name:ident, true) => {
        #[inline]
        fn hash_key(&self) -> String {
            format!(
                "{}({}; [{}]; [{}])",
                stringify!($type_name),
                self.name,
                self.dependencies.hash_key(),
                self.derivative.hash_key(),
            )
        }
    };

    (@impl_hash_key $type_name:ident, false) => {
        #[inline]
        fn hash_key(&self) -> String {
            format!(
                "{}({}; [{}])",
                stringify!($type_name),
                self.name,
                self.derivative.hash_key(),
            )
        }
    };

    (@impl_match_for_find_all $self:ident, $op:ident, true) => {
        $self.name == $op.name
            && $self.dependencies == $op.dependencies
    };

    (@impl_match_for_find_all $self:ident, $op:ident, false) => {
        $self.name == $op.name
    };

    (@impl_eliminate $type_name:ident, true) => { };

    (@impl_eliminate $type_name:ident, false) => {
        #[inline]
        fn eliminate(
            &self,
            parameter: &Arc<dyn Expr>,
            perturbations: &[Arc<Perturbation>],
            min_order: u32,
        ) -> Result<Arc<dyn Expr>, TinnedError>
        {
            if let Some(op) = downcast_from_arc::<$type_name>(parameter) {
                if self.name == op.name {
                    let map = self.derivative.get_map_clone();
                    let order: u32 = perturbations
                        .iter()
                        .map(|p| *map.get(p).unwrap_or(&0))
                        .sum();
                    if order >= min_order && order <= perturbations.len() as u32 {
                        return Ok(ZeroOperator::new());
                    }
                }
            }

            Ok(self.clone_expr())
        }
    };
}

macro_rules! impl_nullary_expr_common_methods {
    ($is_scalar:tt) => {
        impl_expr_common_methods!($is_scalar);

        #[inline]
        fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
            if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                impl_zero_expr!($is_scalar)
            } else {
                Ok(self.clone_expr())
            }
        }

        #[inline]
        fn replace_all(
            &self,
            map: &HashMap<Arc<dyn Expr>, Arc<dyn Expr>>,
        ) -> Result<Arc<dyn Expr>, TinnedError> {
            if let Some((_, value)) = map.iter().find(|(key, _)| self.match_for_replace_all(key)) {
                if self.derivative.is_empty() {
                    Ok(value.clone())
                } else {
                    differentiate_expr(value, &self.derivative)
                }
            } else {
                Ok(self.clone_expr())
            }
        }
    };
}
