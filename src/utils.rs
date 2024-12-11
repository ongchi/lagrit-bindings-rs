use std::path::{Path, PathBuf};

use crate::error::LagritError;

pub struct Pushd {
    current: PathBuf,
}

impl Pushd {
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, LagritError> {
        let current = std::env::current_dir()?;
        let target = path.as_ref().to_path_buf();

        if !target.exists() {
            std::fs::create_dir_all(&target)?;
        }

        std::env::set_current_dir(&target)?;
        Ok(Self { current })
    }

    pub fn pop(self) -> Result<(), LagritError> {
        let _ = std::env::set_current_dir(&self.current);
        Ok(())
    }
}
