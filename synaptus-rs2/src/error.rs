use std::fmt;

/// Errors returned by the library. Numerical failure is never converted into
/// success-shaped output.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// Input violates a documented precondition.
    InvalidInput(String),
    /// A size computation overflowed or exceeded what the algorithm supports.
    SizeOverflow(String),
    /// I/O or parse failure when reading exported datasets.
    Io(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidInput(s) => write!(f, "invalid input: {s}"),
            Error::SizeOverflow(s) => write!(f, "size overflow: {s}"),
            Error::Io(s) => write!(f, "I/O error: {s}"),
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

pub(crate) fn invalid<T>(msg: impl Into<String>) -> Result<T> {
    Err(Error::InvalidInput(msg.into()))
}
