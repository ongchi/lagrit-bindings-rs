use itertools::Itertools;
use std::{
    ffi::CString,
    fs::File,
    io::{BufReader, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering},
        Arc, OnceLock, RwLock,
    },
};

use lagrit_sys::{fc_fclose, fc_fflush_and_sync};

use super::{InitMode, MeshObject};
use crate::{
    error::LagritError,
    ffi::{cmo_get_name, dotask, initlagrit, mmfindbk_string, mmrelprt},
    utils::Pushd,
};

pub(crate) static LAGRIT: OnceLock<Arc<LaGriT>> = OnceLock::new();

pub struct LaGriT {
    is_initialized: AtomicBool,
    log_file: RwLock<String>,
    batch_file: RwLock<String>,
    workdir: RwLock<PathBuf>,
    cmd_msg: RwLock<String>,
    msg_last_pos: AtomicUsize,
}

impl LaGriT {
    pub fn new(
        mode: InitMode,
        log_file: Option<&str>,
        batch_file: Option<&str>,
        workdir: Option<&str>,
    ) -> Result<Arc<Self>, LagritError> {
        let lg = LAGRIT.get_or_init(|| {
            Arc::new(Self {
                is_initialized: AtomicBool::new(false),
                log_file: RwLock::new(String::with_capacity(64)),
                batch_file: RwLock::new(String::with_capacity(64)),
                workdir: RwLock::new(PathBuf::new()),
                cmd_msg: RwLock::new(String::new()),
                msg_last_pos: AtomicUsize::new(0),
            })
        });

        if lg.is_initialized.load(Ordering::Relaxed) {
            return Err(LagritError::AlreadyInitialized);
        }

        // Clear last command message
        lg.cmd_msg
            .write()
            .map(|mut s| s.clear())
            .map_err(|_| LagritError::RwLockPoisoned)?;
        lg.msg_last_pos.store(0, Ordering::Relaxed);

        let workdir = if let Some(d) = workdir {
            PathBuf::from(d)
        } else {
            std::env::current_dir()?
        };

        let d = Pushd::new(&workdir)?;

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
        lg.workdir
            .write()
            .map(|mut p| *p = workdir)
            .map_err(|_| LagritError::RwLockPoisoned)?;
        lg.update_message()?;

        d.pop()?;

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
        let workdir = self
            .workdir
            .read()
            .map_err(|_| LagritError::RwLockPoisoned)?
            .clone();
        let d = Pushd::new(&workdir)?;

        dotask(cmd)?;
        self.update_message()?;

        d.pop()?;

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

    pub fn read_mo<P: AsRef<Path>>(
        &self,
        file_path: P,
        name: Option<&str>,
    ) -> Result<Vec<MeshObject>, LagritError> {
        let file_path = file_path.as_ref();
        let file_name = file_path
            .file_name()
            .and_then(|s| s.to_str())
            .ok_or_else(|| LagritError::InvalidPath(file_path.to_string_lossy().to_string()))?;
        let file_ext = file_path
            .extension()
            .and_then(|e| e.to_str())
            .ok_or_else(|| LagritError::InvalidPath(file_path.to_string_lossy().to_string()))?;

        let (is_lg, mo_name) = match file_ext {
            // mo_name is ignored for lagrit file
            "lg" | "lagrit" => (true, "".to_string()),
            _ => (
                false,
                name.map(|n| n.to_string()).unwrap_or_else(|| {
                    file_path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .map(|s| s.to_string())
                        .unwrap_or(self.new_name().unwrap())
                }),
            ),
        };

        let current_mos = self.mo_names()?;

        // Check if file exists
        let file = File::open(file_path)?;
        drop(file);

        let workdir = self
            .workdir
            .read()
            .map_err(|_| LagritError::RwLockPoisoned)?
            .clone();
        let currentdir = std::env::current_dir()?;

        // Create symbolic link if workdir is different from currentdir
        let link_path = if workdir != currentdir {
            let file_path = if file_path.is_absolute() {
                file_path.to_path_buf()
            } else {
                currentdir.join(file_path)
            };
            let link_path = workdir.join(file_name);
            std::os::unix::fs::symlink(file_path, &link_path)?;
            Some(link_path)
        } else {
            None
        };

        if is_lg {
            self.sendcmd(&format!("read/lagrit/{file_name}"))?;
        } else {
            self.sendcmd(&format!("read/{file_name}/{mo_name}"))?;
        }

        // Remove symbolic link
        link_path.map(std::fs::remove_file).transpose()?;

        let new_mos = self.mo_names()?;

        Ok(new_mos
            .into_iter()
            .filter(|n| !current_mos.contains(n))
            .map(|n| MeshObject::new(&n))
            .collect())
    }

    pub fn dump_mo<P: AsRef<Path>>(
        &self,
        file_path: P,
        mo: Option<MeshObject>,
    ) -> Result<(), LagritError> {
        let file_path = file_path.as_ref();
        let file_ext = file_path.extension().and_then(|e| e.to_str());
        let mut file_name = file_path
            .file_name()
            .and_then(|s| s.to_str())
            .map(|s| s.to_string())
            .ok_or_else(|| LagritError::InvalidPath(file_path.to_string_lossy().to_string()))?;

        let file_ext = file_ext.unwrap_or_else(|| {
            file_name.push_str(".lg");
            ".lg"
        });

        if let Some(mo) = mo {
            self.sendcmd(&format!("dump/{file_name}/{}", mo.name()))?;
        } else if file_ext == "lg" || file_ext == "lagrit" {
            self.sendcmd(&format!("dump/{file_name}/-all-"))?;
        } else {
            return Err(LagritError::InvalidArguments(
                "Only one MeshObject is allowed for non-Lagrit file".to_string(),
            ));
        };

        let parent_path = file_path.parent().unwrap_or_else(|| Path::new("."));
        std::fs::create_dir_all(parent_path)?;
        let workdir = self
            .workdir
            .read()
            .map_err(|_| LagritError::RwLockPoisoned)?
            .clone();
        let dump_file = workdir.join(file_name);
        std::fs::rename(&dump_file, file_path)?;

        Ok(())
    }

    fn new_name(&self) -> Result<String, LagritError> {
        let current_mos = self.mo_names()?;
        match current_mos
            .iter()
            .filter_map(|n| {
                if (n.len() > 2) && n.starts_with("mo") {
                    n[2..].parse::<usize>().ok()
                } else {
                    None
                }
            })
            .sorted()
            .last()
        {
            Some(last) => Ok(format!("mo{}", last + 1)),
            None => Ok("mo1".to_string()),
        }
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

        let workdir = self
            .workdir
            .read()
            .map_err(|_| LagritError::RwLockPoisoned)?
            .clone();
        let d = Pushd::new(workdir)?;

        unsafe {
            fc_fclose(log_file.as_ptr());
            fc_fclose(batch_file.as_ptr());
        }

        d.pop()?;

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
