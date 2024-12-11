use std::sync::Arc;

use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::{create_exception, PyErr};

pub use crate::error::LagritError;
use crate::objects::{AttrValue, LaGriT, MeshObject};

create_exception!(pylagrit, PyLagritError, PyException);

impl From<LagritError> for PyErr {
    fn from(e: LagritError) -> Self {
        PyLagritError::new_err(e.to_string())
    }
}

#[pyclass(name = "MeshObject", module = "lagrit_bindings", subclass)]
pub struct PyMeshObject {
    pub mesh_object: Arc<MeshObject>,
}

#[pymethods]
impl PyMeshObject {
    #[new]
    #[pyo3(signature = (name))]
    fn new(name: &str) -> PyResult<Self> {
        Ok(Self {
            mesh_object: Arc::new(MeshObject::new(name)),
        })
    }

    fn name(&self) -> &str {
        self.mesh_object.name()
    }

    fn status(&self) -> PyResult<()> {
        self.mesh_object.status().map_err(|e| e.into())
    }

    fn mesh_type(&self) -> PyResult<(String, i32)> {
        self.mesh_object.mesh_type().map_err(|e| e.into())
    }

    fn attr(&self, attr: &str) -> PyResult<AttrValue> {
        self.mesh_object.attr(attr).map_err(|e| e.into())
    }

    fn attr_list(&self) -> PyResult<Vec<Vec<String>>> {
        self.mesh_object.attr_list().map_err(|e| e.into())
    }

    fn pset_names(&self) -> PyResult<Vec<String>> {
        self.mesh_object.pset_names().map_err(|e| e.into())
    }

    fn eltset_names(&self) -> PyResult<Vec<String>> {
        self.mesh_object.eltset_names().map_err(|e| e.into())
    }
}

#[pyclass(name = "LaGriT", module = "lagrit_bindings", subclass)]
pub struct PyLaGriT {
    pub lagrit: Arc<LaGriT>,
}

#[pymethods]
impl PyLaGriT {
    #[new]
    #[pyo3(signature = (mode, log_file=None, batch_file=None))]
    fn new(mode: &str, log_file: Option<&str>, batch_file: Option<&str>) -> PyResult<Self> {
        Ok(Self {
            lagrit: LaGriT::new(mode.parse()?, log_file, batch_file)?,
        })
    }

    fn sendcmd(&self, cmd: &str) -> PyResult<()> {
        Ok(self.lagrit.sendcmd(cmd)?)
    }

    fn cmdmsg(&self) -> PyResult<String> {
        Ok(self.lagrit.cmdmsg()?)
    }

    fn mo_names(&self) -> PyResult<Vec<String>> {
        Ok(self.lagrit.mo_names()?)
    }

    fn pset_names(&self) -> PyResult<Vec<String>> {
        Ok(self.lagrit.pset_names()?)
    }

    fn fset_names(&self) -> PyResult<Vec<String>> {
        Ok(self.lagrit.fset_names()?)
    }

    fn eltset_names(&self) -> PyResult<Vec<String>> {
        Ok(self.lagrit.eltset_names()?)
    }

    fn cmo(&self) -> PyResult<PyMeshObject> {
        Ok(PyMeshObject {
            mesh_object: Arc::new(self.lagrit.cmo()?),
        })
    }

    fn close(&self) -> PyResult<()> {
        Ok(self.lagrit.close()?)
    }
}

#[pymodule]
fn lagrit_bindings(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMeshObject>()?;
    m.add_class::<PyLaGriT>()?;
    Ok(())
}
