macro_rules! dot_product_argument_op {
    ($self:ident, $bra_operation:expr, $ket_operation:expr, $message:literal) => {{
        let new_bra = $bra_operation.map_err(|e| {
            generic_expression_error(concat!($message, " for bra"), $self, Some(Box::new(e)))
        })?;
        let new_ket = $ket_operation.map_err(|e| {
            generic_expression_error(concat!($message, " for ket"), $self, Some(Box::new(e)))
        })?;

        if &new_bra == &$self.bra && &new_ket == &$self.ket {
            Ok($self.clone_expr())
        } else {
            Self::make_dot_product(new_bra, new_ket, $self.allow_braket_swap)
        }
    }};
}
