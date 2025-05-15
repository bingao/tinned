macro_rules! impl_adjoint_map_operation {
    ($self:ident, $operation:expr, $message:expr) => {{
        let new_target = ($operation)(&$self.target).map_err(|e| {
            generic_expression_error(concat!($message, " for target"), $self, Some(Box::new(e)))
        })?;

        let mut new_adj_map = &new_target != &$self.target;

        let mut new_adj_chain = Vec::with_capacity($self.adjoint_chain.len());

        for x in &$self.adjoint_chain {
            let new_x = ($operation)(x).map_err(|e| {
                generic_expression_error(concat!($message, " for chain"), $self, Some(Box::new(e)))
            })?;
            if is_zero_expr(&new_x, None) {
                return Ok(ZeroOperator::new());
            } else {
                if !new_adj_map {
                    new_adj_map = &new_x != x;
                }
                new_adj_chain.push(new_x);
            }
        }

        if new_adj_map {
            Ok(Self::new(
                new_adj_chain,
                $self.chain_commutative,
                new_target,
                Some($self.left_action),
            ))
        } else {
            Ok($self.clone_expr())
        }
    }};
}
