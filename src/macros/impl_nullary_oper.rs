macro_rules! impl_nullary_oper_type {
    ($type_name:ident, $builder_name:ident, $has_deps:tt, $is_scalar:tt) => {
        impl_nullary_oper_type!(@def_oper_struct $type_name, $has_deps);

        impl $type_name {
            impl_nullary_oper_type!(@impl_oper_methods $builder_name, $has_deps);

            #[inline]
            pub fn name(&self) -> &str {
                &self.name
            }

            #[inline]
            pub fn derivative(&self) -> &PertMultichain {
                &self.derivative
            }
        }

        impl_nullary_oper_type!(@def_builder_struct $builder_name, $has_deps);

        impl $builder_name {
            #[inline]
            pub fn derivative(mut self, deriv: PertMultichain) -> Self {
                self.derivative = deriv;
                self
            }

            impl_nullary_oper_type!(@impl_builder_methods $type_name, $has_deps, $is_scalar);
        }
    };

    (@def_oper_struct $type_name:ident, true) => {
        #[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
        pub struct $type_name {
            name: String,
            dependencies: PertMultichain,
            derivative: PertMultichain,
        }
    };

    (@def_oper_struct $type_name:ident, false) => {
        #[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
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
            if is_sub_multichain(&self.derivative, &self.dependencies) {
                Ok(intern(Arc::new($type_name {
                    name: self.name,
                    dependencies: self.dependencies,
                    derivative: self.derivative,
                })))
            } else {
                Ok(0.into())
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
            if is_sub_multichain(&self.derivative, &self.dependencies) {
                Ok(intern(Arc::new($type_name {
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
        compile_error!("impl_nullary_oper_type!(...) does not support has_deps = false and is_scalar = true");
    };

    (@impl_builder_methods $type_name:ident, false, false) => {
        #[inline]
        pub fn build(self) -> Result<Arc<dyn Expr>, TinnedError> {
            Ok(intern(Arc::new($type_name {
                name: self.name,
                derivative: self.derivative,
            })))
        }
    };
}

macro_rules! impl_nullary_oper_traits {
    ($type_name:ident, $has_deps:tt, $is_scalar:tt) => {
        impl Expr for $type_name {
            #[inline]
            fn as_any(&self) -> &dyn Any {
                self
            }

            impl_nullary_oper_traits!(@impl_hash_key $type_name, $has_deps);

            #[inline]
            fn is_scalar(&self) -> bool {
                $is_scalar
            }

            #[inline]
            fn eq_expr(&self, other: &dyn Expr) -> bool {
                if let Some(op) = downcast_expr::<$type_name>(other) {
                    impl_nullary_oper_traits!(@impl_eq_expr $has_deps, self, op)
                } else {
                    false
                }
            }

            #[inline]
            fn differentiate(&self, s: &Perturbation) -> Result<Arc<dyn Expr>, TinnedError>
            {
                let mut new_deriv = self.derivative.clone();
                *new_deriv.entry(s.clone()).or_insert(0) += 1;

                self.builder_from(new_deriv).build()
            }
        }

        impl Display for $type_name {
            fn fmt(&self, f: &mut Formatter) -> FmtResult {
                if self.derivative.is_empty() {
                    write!(f, "{}", self.name)
                } else {
                    write!(
                        f,
                        "{}^({})",
                        self.name,
                        pert_multichain_display(&self.derivative),
                    )
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
                pert_multichain_hash_key(&self.dependencies),
                pert_multichain_hash_key(&self.derivative),
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
                pert_multichain_hash_key(&self.derivative),
            )
        }
    };

    (@impl_eq_expr true, $self_:ident, $op:ident) => {
        $self_.name == $op.name
            && $self_.dependencies == $op.dependencies
            && $self_.derivative == $op.derivative
    };

    (@impl_eq_expr false, $self_:ident, $op:ident) => {
        $self_.name == $op.name
            && $self_.derivative == $op.derivative
    };
}
