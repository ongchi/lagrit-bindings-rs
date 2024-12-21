use std::sync::Arc;

use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::{PyFloat, PyInt};
use pyo3::{create_exception, PyErr};

pub use crate::error::LagritError;
use crate::objects::{AttrInfo, AttrValue, CmdWithInput, CmdWithOutput, LaGriT, MeshObject};

create_exception!(pylagrit, PyLagritError, PyException);

impl From<LagritError> for PyErr {
    fn from(e: LagritError) -> Self {
        PyLagritError::new_err(e.to_string())
    }
}

impl<'py> FromPyObject<'py> for AttrValue {
    fn extract_bound(obj: &Bound<'py, PyAny>) -> PyResult<Self> {
        if obj.downcast::<PyInt>().is_ok() {
            Ok(AttrValue::Int(obj.extract()?))
        } else if obj.downcast::<PyFloat>().is_ok() {
            Ok(AttrValue::Real(obj.extract()?))
        } else {
            match obj.extract::<Vec<i64>>() {
                Ok(v) => Ok(AttrValue::Vint(v)),
                Err(_) => Ok(AttrValue::Vdouble(obj.extract()?)),
            }
        }
    }
}

impl<'py> FromPyObject<'py> for MeshObject {
    fn extract_bound(obj: &Bound<'py, PyAny>) -> PyResult<Self> {
        match obj.downcast::<PyMeshObject>() {
            Ok(mo) => Ok(MeshObject::new(mo.extract::<PyMeshObject>()?.name())),
            Err(_) => Err(PyErr::new::<PyLagritError, _>("Invalid MeshObject")),
        }
    }
}

#[derive(Clone)]
#[pyclass(name = "MeshObject", module = "lagrit_bindings", subclass)]
pub struct PyMeshObject {
    pub mesh_object: Arc<MeshObject>,
}

impl From<MeshObject> for PyMeshObject {
    fn from(mo: MeshObject) -> Self {
        Self {
            mesh_object: Arc::new(mo),
        }
    }
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

    fn set_attr(&self, attr: &str, value: AttrValue) -> PyResult<()> {
        self.mesh_object.set_attr(attr, value).map_err(|e| e.into())
    }

    fn attr_list(&self) -> PyResult<Vec<AttrInfo>> {
        self.mesh_object.attr_list().map_err(|e| e.into())
    }

    fn pset_names(&self) -> PyResult<Vec<String>> {
        self.mesh_object.pset_names().map_err(|e| e.into())
    }

    fn eltset_names(&self) -> PyResult<Vec<String>> {
        self.mesh_object.eltset_names().map_err(|e| e.into())
    }

    fn dump(&self, file_path: &str) -> PyResult<()> {
        self.mesh_object.dump(file_path).map_err(|e| e.into())
    }
}

#[pyclass(name = "LgCmdWithInput", module = "lagrit_bindings", subclass)]
pub struct LgCmdWithInput(Arc<CmdWithInput>);

#[pymethods]
impl LgCmdWithInput {
    fn sendcmd(&self, cmd: &str) -> PyResult<()> {
        Ok(self.0.sendcmd(cmd)?)
    }
}

#[pyclass(name = "LgCmdWithOutput", module = "lagrit_bindings", subclass)]
pub struct LgCmdWithOutput(Arc<CmdWithOutput>);

#[pymethods]
impl LgCmdWithOutput {
    fn sendcmd(&self, cmd: &str) -> PyResult<()> {
        Ok(self.0.sendcmd(cmd)?)
    }
}

#[pyclass(name = "LaGriT", module = "lagrit_bindings", subclass)]
pub struct PyLaGriT {
    pub lagrit: Arc<LaGriT>,
}

#[pymethods]
impl PyLaGriT {
    #[new]
    #[pyo3(signature = (mode, log_file=None, batch_file=None, workdir=None))]
    fn new(
        mode: &str,
        log_file: Option<&str>,
        batch_file: Option<&str>,
        workdir: Option<&str>,
    ) -> PyResult<Self> {
        Ok(Self {
            lagrit: LaGriT::new(mode.parse()?, log_file, batch_file, workdir)?,
        })
    }

    fn sendcmd(&self, cmd: &str) -> PyResult<()> {
        Ok(self.lagrit.sendcmd(cmd)?)
    }

    fn with_input(&self, input: &str) -> PyResult<LgCmdWithInput> {
        Ok(LgCmdWithInput(Arc::new(self.lagrit.with_input(input))))
    }

    fn with_output(&self, output: &str) -> PyResult<LgCmdWithOutput> {
        Ok(LgCmdWithOutput(Arc::new(self.lagrit.with_output(output))))
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

    #[pyo3(signature = (name=None))]
    fn get_mo(&self, name: Option<&str>) -> PyResult<PyMeshObject> {
        Ok(PyMeshObject {
            mesh_object: Arc::new(self.lagrit.get_mo(name)?),
        })
    }

    #[pyo3(signature = (file_path, name=None))]
    fn read_mo(&self, file_path: &str, name: Option<&str>) -> PyResult<Vec<PyMeshObject>> {
        Ok(self
            .lagrit
            .read_mo(file_path, name)?
            .into_iter()
            .map(|mo| mo.into())
            .collect())
    }

    #[pyo3(signature = (file_path, mo=None))]
    fn dump_mo(&self, file_path: &str, mo: Option<MeshObject>) -> PyResult<()> {
        Ok(self.lagrit.dump_mo(file_path, mo)?)
    }

    fn close(&self) -> PyResult<()> {
        Ok(self.lagrit.close()?)
    }
}

#[pymodule]
fn lagrit_bindings(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMeshObject>()?;
    m.add_class::<LgCmdWithInput>()?;
    m.add_class::<LgCmdWithOutput>()?;
    m.add_class::<PyLaGriT>()?;
    Ok(())
}
