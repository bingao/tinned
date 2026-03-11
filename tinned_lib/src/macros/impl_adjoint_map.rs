macro_rules! impl_adjoint_map_operation {
    ($self:ident, $operation:expr, $message:expr) => {{
        let new_target = ($operation)(&$self.target).map_err(|e| {
            $crate::public::generic_expression_error(
                concat!($message, " for target"),
                $self,
                Some(::std::boxed::Box::new(e)),
            )
        })?;
        if $crate::public::is_zero_expr(&new_target, None) {
            return Ok($crate::expressions::ZeroOperator::new());
        }

        let mut new_ad_map = &new_target != &$self.target;

        let mut new_generators = ::std::vec::Vec::with_capacity($self.generators.len());

        for generator in &$self.generators {
            let new_generator = ($operation)(generator).map_err(|e| {
                $crate::public::generic_expression_error(
                    concat!($message, " for generators"),
                    $self,
                    Some(::std::boxed::Box::new(e)),
                )
            })?;
            if $crate::public::is_zero_expr(&new_generator, None) {
                return Ok($crate::expressions::ZeroOperator::new());
            } else {
                if !new_ad_map {
                    new_ad_map = &new_generator != generator;
                }
                new_generators.push(new_generator);
            }
        }

        if new_ad_map {
            Self::new(new_generators, new_target, Some($self.left_action))
        } else {
            Ok($self.clone_expr())
        }
    }};
}
