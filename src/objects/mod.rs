pub mod lg;
pub mod mo;

pub use lg::{CmdWithInput, CmdWithOutput, LaGriT};
pub use mo::MeshObject;

use crate::error::LagritError;

#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum InitMode {
    Slient,
    Noisy,
}

impl std::str::FromStr for InitMode {
    type Err = LagritError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "slient" => Ok(Self::Slient),
            "noisy" => Ok(Self::Noisy),
            _ => Err(LagritError::InvalidInitMode),
        }
    }
}

#[derive(Debug)]
#[cfg_attr(feature = "pyo3", derive(pyo3::IntoPyObject))]
pub enum AttrValue {
    Int(i64),
    Vint(Vec<i64>),
    Real(f64),
    Vdouble(Vec<f64>),
    Character(String),
    Vchar(Vec<String>),
}

#[derive(Debug)]
#[cfg_attr(feature = "pyo3", derive(pyo3::IntoPyObject))]
pub struct AttrInfo {
    pub name: String,
    pub type_: String,
    pub rank: String,
    pub length: String,
    pub interpolation: String,
    pub persistence: String,
    pub ioflag: String,
}
