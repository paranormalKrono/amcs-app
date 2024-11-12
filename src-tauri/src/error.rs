use std::sync::PoisonError;

// Create a custom Error that we can return in Results
#[derive(Debug, thiserror::Error)]
pub enum Error {
    // Implement std::io::Error for our Error enum
    #[error(transparent)]
    Io(#[from] std::io::Error),
    // Add a PoisonError, but we implement it manually later
    #[error("the mutex was poisoned")]
    Poison(String),
    #[error("problem with writing in csv file")]
    Scv(String),
}
// Implement Serialize for the error
impl serde::Serialize for Error {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer,
    {
        serializer.serialize_str(self.to_string().as_ref())
    }
}
// Implement From<PoisonError> for Error to convert it to something we have set up serialization for
impl<T> From<PoisonError<T>> for Error {
    fn from(err: PoisonError<T>) -> Self {
        // We "just" convert the error to a string here
        Error::Poison(err.to_string())
    }
}

impl From<csv::Error> for Error {
    fn from(err: csv::Error) -> Self {
        // We "just" convert the error to a string here
        Error::Scv(err.to_string())
    }
}
