use std::{
    ffi::CString,
    fs::File,
    io::{BufReader, Read, Seek, SeekFrom},
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering},
        Arc, OnceLock, RwLock,
    },
};

use lagrit_sys::{fc_fclose, fc_fflush_and_sync};

use crate::{
    error::LagritError,
    ffi::{
        cmo_attlist, cmo_get_info, cmo_get_mesh_type, cmo_get_name, dotask, initlagrit,
        mmfindbk_string, mmrelprt,
    },
};

static LAGRIT: OnceLock<Arc<LaGriT>> = OnceLock::new();

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

pub struct LaGriT {
    is_initialized: AtomicBool,
    log_file: RwLock<String>,
    batch_file: RwLock<String>,
    cmd_msg: RwLock<String>,
    msg_last_pos: AtomicUsize,
}

pub struct MeshObject {
    name: String,
}

pub enum InitMode {
    Slient,
    Noisy,
}

impl LaGriT {
    pub fn new(
        mode: InitMode,
        log_file: Option<&str>,
        batch_file: Option<&str>,
    ) -> Result<Arc<Self>, LagritError> {
        let lg = LAGRIT.get_or_init(|| {
            Arc::new(Self {
                is_initialized: AtomicBool::new(false),
                log_file: RwLock::new(String::with_capacity(64)),
                batch_file: RwLock::new(String::with_capacity(64)),
                cmd_msg: RwLock::new(String::new()),
                msg_last_pos: AtomicUsize::new(0),
            })
        });

        if lg.is_initialized.load(Ordering::Relaxed) {
            return Err(LagritError::AlreadyInitialized);
        }

        lg.cmd_msg
            .write()
            .map(|mut s| s.clear())
            .map_err(|_| LagritError::RwLockPoisoned)?;
        lg.msg_last_pos.store(0, Ordering::Relaxed);

        let log_file = log_file.unwrap_or("lagrit.log");
        let batch_file = batch_file.unwrap_or("lagrit.out");
        initlagrit(
            match mode {
                InitMode::Slient => 0,
                InitMode::Noisy => 1,
            },
            log_file,
            batch_file,
        )?;

        lg.is_initialized.store(true, Ordering::Relaxed);
        lg.log_file
            .write()
            .map(|mut s| {
                s.clear();
                s.push_str(log_file);
            })
            .map_err(|_| LagritError::RwLockPoisoned)?;
        lg.batch_file
            .write()
            .map(|mut s| {
                s.clear();
                s.push_str(batch_file);
            })
            .map_err(|_| LagritError::RwLockPoisoned)?;
        lg.update_message()?;

        Ok(lg.clone())
    }

    fn update_message(&self) -> Result<(), LagritError> {
        let log_file = self
            .log_file
            .read()
            .map_err(|_| LagritError::RwLockPoisoned)?
            .clone();
        let batch_file = self
            .batch_file
            .read()
            .map_err(|_| LagritError::RwLockPoisoned)?
            .clone();

        unsafe {
            fc_fflush_and_sync(log_file.as_ptr() as *mut i8);
            fc_fflush_and_sync(batch_file.as_ptr() as *mut i8);
        }

        let msg_file = File::open(batch_file)?;
        let mut reader = BufReader::new(msg_file);
        reader.seek(SeekFrom::Start(
            self.msg_last_pos.load(Ordering::Relaxed) as u64
        ))?;

        let mut buf = String::new();
        let pos = reader.read_to_string(&mut buf)?;

        self.cmd_msg
            .write()
            .map(|mut s| {
                s.clear();
                s.push_str(&buf);
            })
            .map_err(|_| LagritError::RwLockPoisoned)?;

        self.msg_last_pos.fetch_add(pos, Ordering::Relaxed);

        Ok(())
    }

    pub fn sendcmd(&self, cmd: &str) -> Result<(), LagritError> {
        dotask(cmd)?;
        self.update_message()?;

        Ok(())
    }

    pub fn cmdmsg(&self) -> Result<String, LagritError> {
        Ok(self
            .cmd_msg
            .read()
            .map_err(|_| LagritError::RwLockPoisoned)?
            .clone())
    }

    pub fn mo_names(&self) -> Result<Vec<String>, LagritError> {
        Ok(mmfindbk_string("cmo_names", "define_cmo_lg")?
            .into_iter()
            .filter(|c| !c.is_empty() && c != "-default-")
            .collect())
    }

    pub fn pset_names(&self) -> Result<Vec<String>, LagritError> {
        let mut names = vec![];
        for mo in self.mo_names()? {
            names.extend(MeshObject::new(&mo).pset_names()?);
        }
        Ok(names)
    }

    pub fn fset_names(&self) -> Result<Vec<String>, LagritError> {
        let mut names = vec![];
        for mo in self.mo_names()? {
            names.extend(MeshObject::new(&mo).fset_names()?);
        }
        Ok(names)
    }

    pub fn eltset_names(&self) -> Result<Vec<String>, LagritError> {
        let mut names = vec![];
        for mo in self.mo_names()? {
            names.extend(MeshObject::new(&mo).eltset_names()?);
        }
        Ok(names)
    }

    // Get current MeshObject
    pub fn cmo(&self) -> Result<MeshObject, LagritError> {
        Ok(MeshObject::new(&cmo_get_name()?))
    }

    // Close log and batch files, then release memory
    pub fn close(&self) -> Result<(), LagritError> {
        if !self.is_initialized.load(Ordering::Relaxed) {
            return Err(LagritError::NotInitialized);
        }

        let log_file = CString::new(
            self.log_file
                .read()
                .map_err(|_| LagritError::RwLockPoisoned)?
                .clone(),
        )?;
        let batch_file = CString::new(
            self.batch_file
                .read()
                .map_err(|_| LagritError::RwLockPoisoned)?
                .clone(),
        )?;

        unsafe {
            fc_fclose(log_file.as_ptr());
            fc_fclose(batch_file.as_ptr());
        }

        self.mo_names()
            .unwrap_or_else(|_| vec![])
            .iter()
            .map(|n| n.as_str())
            .chain([
                "global_lg",
                "geom_lg",
                "default_cmo_lg",
                "initlagrit",
                "define_cmo_lg",
            ])
            .for_each(|part| mmrelprt(part).unwrap());

        self.is_initialized.store(false, Ordering::Relaxed);

        Ok(())
    }
}

impl MeshObject {
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn status(&self) -> Result<(), LagritError> {
        dotask(&format!("cmo/status/{}", self.name))
    }

    pub fn mesh_type(&self) -> Result<(String, i32), LagritError> {
        cmo_get_mesh_type(&self.name)
    }

    pub fn attr(&self, attr: &str) -> Result<AttrValue, LagritError> {
        cmo_get_info(attr, &self.name)
    }

    pub fn attr_list(&self) -> Result<Vec<Vec<String>>, LagritError> {
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
}
