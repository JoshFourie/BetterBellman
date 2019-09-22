use std::error::Error;
use std::fmt;
use std::io;
use std::option;

pub type Result<T> = std::result::Result<T, SynthesisError>;

/// This is an error that could occur during circuit synthesis contexts,
/// such as CRS generation, proving or verification.
#[derive(Debug)]
pub enum SynthesisError {
    /// During synthesis, we lacked knowledge of a variable assignment.
    AssignmentMissing,
    /// During synthesis, we divided by zero.
    DivisionByZero,
    /// During synthesis, we constructed an unsatisfiable constraint system.
    Unsatisfiable,
    /// During synthesis, our polynomials ended up being too high of degree
    PolynomialDegreeTooLarge,
    /// During proof generation, we encountered an identity in the CRS
    UnexpectedIdentity,
    /// During proof generation, we encountered an I/O error with the CRS
    IoError(io::Error),
    /// During verification, our verifying key was malformed.
    MalformedVerifyingKey,
    /// During CRS generation, we observed a malformed wire size.
    MalformedWireSize,
    /// During CRS generation, we observed an unconstrained auxiliary variable
    UnconstrainedVariable,
    /// During synthesis, we called an operation on a None.
    Null
}

impl From<option::NoneError> for SynthesisError {
    fn from(_: option::NoneError) -> Self {
        SynthesisError::Null
    }
}

impl From<io::Error> for SynthesisError {
    fn from(e: io::Error) -> SynthesisError {
        SynthesisError::IoError(e)
    }
}

impl Error for SynthesisError {
    fn description(&self) -> &str {
        match *self {
            SynthesisError::AssignmentMissing => {
                "an assignment for a variable could not be computed"
            }
            SynthesisError::DivisionByZero => "division by zero",
            SynthesisError::Unsatisfiable => "unsatisfiable constraint system",
            SynthesisError::PolynomialDegreeTooLarge => "polynomial degree is too large",
            SynthesisError::UnexpectedIdentity => "encountered an identity element in the CRS",
            SynthesisError::IoError(_) => "encountered an I/O error",
            SynthesisError::MalformedVerifyingKey => "malformed verifying key",
            SynthesisError::MalformedWireSize => "malformed wire size",
            SynthesisError::UnconstrainedVariable => "auxiliary variable was unconstrained",
            SynthesisError::Null => "encountered an operation on a None"
        }
    }
}

impl fmt::Display for SynthesisError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let &SynthesisError::IoError(ref e) = self {
            write!(f, "I/O error: ")?;
            e.fmt(f)
        } else {
            write!(f, "{}", self.description())
        }
    }
}
