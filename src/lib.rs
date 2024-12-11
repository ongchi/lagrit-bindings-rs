pub mod error;
pub mod ffi;
pub mod objects;
pub(crate) mod utils;

pub use error::LagritError;
pub use objects::{lg::LaGriT, mo::MeshObject, InitMode};

#[cfg(feature = "pyo3")]
pub mod py;
