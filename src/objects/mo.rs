use super::lg::LAGRIT;
use super::{AttrInfo, AttrValue};
use crate::error::LagritError;
use crate::ffi::{cmo_attlist, cmo_get_info, cmo_get_mesh_type, cmo_set_info, mmfindbk_string};

#[derive(Debug, Clone)]
pub struct MeshObject {
    name: String,
}

impl MeshObject {
    pub(crate) fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn status(&self) -> Result<(), LagritError> {
        let lg = LAGRIT.get().ok_or(LagritError::NotInitialized)?;
        if lg.is_initialized() {
            lg.sendcmd(&format!("cmo/status/{}", self.name))
        } else {
            Err(LagritError::NotInitialized)
        }
    }

    pub fn mesh_type(&self) -> Result<(String, i32), LagritError> {
        cmo_get_mesh_type(&self.name)
    }

    pub fn attr(&self, attr: &str) -> Result<AttrValue, LagritError> {
        cmo_get_info(attr, &self.name)
    }

    pub fn set_attr(&self, attr: &str, value: AttrValue) -> Result<(), LagritError> {
        cmo_set_info(attr, &self.name, value)
    }

    pub fn attr_list(&self) -> Result<Vec<AttrInfo>, LagritError> {
        cmo_attlist(&self.name)
    }

    pub fn pset_names(&self) -> Result<Vec<String>, LagritError> {
        Ok(mmfindbk_string("psetnames", &self.name)?
            .into_iter()
            .filter(|c| !c.is_empty())
            .collect())
    }

    pub fn fset_names(&self) -> Result<Vec<String>, LagritError> {
        Ok(mmfindbk_string("fsetnames", &self.name)?
            .into_iter()
            .filter(|c| !c.is_empty())
            .collect())
    }

    pub fn eltset_names(&self) -> Result<Vec<String>, LagritError> {
        Ok(mmfindbk_string("eltsetnames", &self.name)?
            .into_iter()
            .filter(|c| !c.is_empty())
            .collect())
    }

    pub fn dump(&self, file_path: &str) -> Result<(), LagritError> {
        let lg = LAGRIT.get().ok_or(LagritError::NotInitialized)?;
        lg.dump_mo(file_path, Some(self.clone()))
    }
}
