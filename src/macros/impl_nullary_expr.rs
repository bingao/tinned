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
        #[typetag::serde]
        impl Expr for $type_name {
            #[inline]
            fn as_any(&self) -> &dyn std::any::Any {
                self
            }

            impl_nullary_expr_traits!(@impl_hash_key $type_name, $has_deps);

            #[inline]
            fn is_scalar(&self) -> bool {
                $is_scalar
            }

            #[inline]
            fn clone_expr(&self) -> Arc<dyn Expr> {
                Arc::new(self.clone())
            }

            #[inline]
            fn eq_expr(&self, other: &dyn Expr) -> bool {
                if let Some(op) = downcast_from_ref::<$type_name>(other) {
                    impl_nullary_expr_traits!(@impl_eq_expr self, op, $has_deps)
                } else {
                    false
                }
            }

            #[inline]
            fn eq_shallow(&self, other: &Arc<dyn Expr>) -> bool {
                if let Some(op) = downcast_from_arc::<$type_name>(other) {
                    impl_nullary_expr_traits!(@impl_eq_shallow self, op, $has_deps)
                } else {
                    false
                }
            }

            #[inline]
            fn fmt_expr(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "{self}")
            }

            #[inline]
            fn total_order(&self) -> u32 { self.derivative.total_order() }

            #[inline]
            fn differentiate(&self, s: &Arc<Perturbation>) -> Result<Arc<dyn Expr>, TinnedError>
            {
                let new_deriv = self.derivative.clone_with_insert(s);

                self.builder_from(new_deriv).build()
            }

            impl_nullary_expr_traits!(@impl_eliminate $type_name, $has_deps);

            #[inline]
            fn remove(&self, set: &HashSet<Arc<dyn Expr>>) -> Result<Arc<dyn Expr>, TinnedError> {
                if set.iter().any(|expr| self.eq_expr(expr.as_ref())) {
                    impl_nullary_expr_traits!(@build_zero_expr $is_scalar)
                } else {
                    Ok(self.clone_expr())
                }
            }
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

    (@impl_eq_expr $self:ident, $op:ident, true) => {
        $self.name == $op.name
            && $self.dependencies == $op.dependencies
            && $self.derivative == $op.derivative
    };

    (@impl_eq_expr $self:ident, $op:ident, false) => {
        $self.name == $op.name
            && $self.derivative == $op.derivative
    };

    (@impl_eq_shallow $self:ident, $op:ident, true) => {
        $self.name == $op.name
            && $self.dependencies == $op.dependencies
    };

    (@impl_eq_shallow $self:ident, $op:ident, false) => {
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

    (@build_zero_expr true) => { Ok(Number::zero()) };

    (@build_zero_expr false) => { Ok(ZeroOperator::new()) };
}
