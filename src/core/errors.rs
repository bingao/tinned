/// Custom error type for the Tinned symbolic algebra library
#[derive(Debug, Clone)]
pub enum TinnedError {
    /// Encounter an invalid expression
    InvalidExpression {
        message: &'static str,
        expression: String,
    },

    /// Division by zero
    DivisionByZero,

    /// Error that a code path should never be hit under correct logic
    Unreachable {
        message: &'static str,
        expression: String,
    },

    /// Unimplemented behavior
    NotYetImplemented(String),

    /// Generic message error
    Message(String),
}

impl std::fmt::Display for TinnedError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            TinnedError::InvalidExpression {
                message,
                expression,
            } => {
                write!(f, "Invalid expression: {} encountered: {}", expression, message)
            },
            TinnedError::DivisionByZero => f.write_str("Division by zero encountered"),
            TinnedError::Unreachable {
                message,
                expression,
            } => {
                write!(f, "Unreachable error: {} in expression: {}", message, expression)
            },
            TinnedError::NotYetImplemented(feature) => {
                write!(f, "Feature not implemented: {}", feature)
            },
            TinnedError::Message(msg) => write!(f, "{}", msg),
        }
    }
}

impl std::error::Error for TinnedError {}
