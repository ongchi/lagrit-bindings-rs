pub mod error;
pub mod ffi;
pub mod obj;

pub use error::LagritError;
pub use obj::{InitMode, LaGriT, MeshObject};

#[cfg(feature = "pyo3")]
pub mod py;
