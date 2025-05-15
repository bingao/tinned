macro_rules! impl_adjoint_map_operation {
    ($self:ident, $operation:expr, $message:expr) => {{
        let new_target = ($operation)(&$self.target).map_err(|e| {
            generic_expression_error(concat!($message, " for target"), $self, Some(Box::new(e)))
        })?;

        let mut new_ad_map = &new_target != &$self.target;

        let mut new_generators = Vec::with_capacity($self.generators.len());

        for generator in &$self.generators {
            let new_generator = ($operation)(generator).map_err(|e| {
                generic_expression_error(
                    concat!($message, " for generators"),
                    $self,
                    Some(Box::new(e)),
                )
            })?;
            if is_zero_expr(&new_generator, None) {
                return Ok(ZeroOperator::new());
            } else {
                if !new_ad_map {
                    new_ad_map = &new_generator != generator;
                }
                new_generators.push(new_generator);
            }
        }

        if new_ad_map {
            Ok(Self::new(
                new_generators,
                $self.generator_commutative,
                new_target,
                Some($self.left_action),
            ))
        } else {
            Ok($self.clone_expr())
        }
    }};
}
