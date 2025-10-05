use std::error::Error;

#[derive(Debug, thiserror::Error)]
pub enum TinnedError {
    #[error("Expression error: {message} in expression: {expression}")]
    ExpressionError {
        message: &'static str,
        expression: String,

        #[source]
        source: Option<Box<dyn Error + Send + Sync>>,
    },

    #[error("Perturbation error: {message} in perturbation: {perturbation}")]
    PerturbationError {
        message: &'static str,
        perturbation: String,

        #[source]
        source: Option<Box<dyn Error + Send + Sync>>,
    },

    #[error("Unreachable code: {message} in expression: {expression}")]
    Unreachable {
        message: &'static str,
        expression: String,

        #[source]
        source: Option<Box<dyn Error + Send + Sync>>,
    },

    #[error("{message}")]
    GenericError {
        message: String,

        #[source]
        source: Option<Box<dyn Error + Send + Sync>>,
    },
}
