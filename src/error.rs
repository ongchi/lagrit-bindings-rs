#[derive(thiserror::Error, Debug)]
pub enum LagritError {
    #[error("Invalid initialization mode")]
    InvalidInitMode,
    #[error("LaGriT is already initialized")]
    AlreadyInitialized,
    #[error("LaGriT is not initialized")]
    NotInitialized,
    #[error("RwLock is poisoned")]
    RwLockPoisoned,
    #[error("Unknown error")]
    Unknown,
    #[error("Invalid arguments: {0}")]
    InvalidArguments(String),
    #[error("Zero length result")]
    ZeroLengthResult,
    #[error("Mesh object not found")]
    MeshObjectNotFound,
    #[error("Attribute not found")]
    AttributeNotFound,
    #[error("Unexpected attribute length: {0}")]
    UnexpectLengthValue(String),
    #[error("Unsupported data type")]
    UnsupportedDataType,
    #[error("Invalid value length")]
    InvalidValueLength,
    #[error("Invalid value type")]
    InvalidValueType,
    #[error("Invalid path: {0}")]
    InvalidPath(String),
    #[error("cmd: {0}\n{1}")]
    ErrorWithMessage(String, String),
    #[error("Error code: {0}")]
    Code(i32),
    #[error(transparent)]
    NullError(#[from] std::ffi::NulError),
    #[error(transparent)]
    StdIoError(#[from] std::io::Error),
}

impl From<i32> for LagritError {
    fn from(err_code: i32) -> Self {
        match err_code {
            1 => LagritError::Unknown,
            2 => LagritError::InvalidArguments("".to_string()),
            _ => LagritError::Code(err_code),
        }
    }
}
