use std::ffi::CString;

const STR_BUF_LEN: usize = 32;

use crate::error::LagritError;
use crate::objects::{AttrInfo, AttrValue};
use lagrit_sys::{
    fc_attr_len, fc_cmo_get_index, fc_cmo_get_mesh_type, fc_cmo_get_name, fc_dotask, fc_initlagrit,
    fc_mmfindbk, fc_mmrelprt, fc_set_fattr, fc_set_iattr,
};

pub fn initlagrit(mode: i32, log_file: &str, batch_file: &str) -> Result<(), LagritError> {
    let log_file = CString::new(log_file)?;
    let batch_file = CString::new(batch_file)?;

    unsafe {
        fc_initlagrit(
            (&mode) as *const i32 as *mut i32,
            log_file.as_ptr(),
            batch_file.as_ptr(),
        );
    };

    Ok(())
}

pub fn dotask(cmd: &str) -> Result<(), LagritError> {
    let cmd = CString::new(cmd)?;
    let mut status = 0;
    unsafe { fc_dotask(cmd.as_ptr(), &mut status) };

    if status != 0 {
        return Err(status.into());
    }

    Ok(())
}

pub fn mmfindbk(
    byte_len: usize,
    cell_len: usize,
    aname: &str,
    pname: &str,
) -> Result<Vec<u8>, LagritError> {
    let aname = CString::new(aname)?;
    let pname = CString::new(pname)?;

    let mut arr_len = 0;
    let mut status = 0;

    unsafe { fc_attr_len(aname.as_ptr(), pname.as_ptr(), &mut arr_len, &mut status) };

    if status != 0 {
        return Err(status.into());
    }
    if arr_len == 0 {
        return Err(LagritError::ZeroLengthResult);
    }

    let mut buffer: Vec<u8> = Vec::with_capacity(byte_len * cell_len * (arr_len as usize));

    unsafe {
        fc_mmfindbk(
            &byte_len as *const usize as *mut usize,
            &cell_len as *const usize as *mut usize,
            aname.as_ptr(),
            pname.as_ptr(),
            buffer.as_ptr() as *mut i8,
            &mut arr_len,
            &mut status,
        );
        buffer.set_len(byte_len * cell_len * (arr_len as usize));
    };

    if status != 0 {
        Err(status.into())
    } else {
        Ok(buffer)
    }
}

pub fn mmfindbk_i64(aname: &str, pname: &str) -> Result<Vec<i64>, LagritError> {
    let buffer = mmfindbk(8, 1, aname, pname)?;
    let buffer_i64 =
        unsafe { std::slice::from_raw_parts(buffer.as_ptr() as *const i64, buffer.len() / 8) };
    Ok(buffer_i64.to_vec())
}

pub fn mmfindbk_f64(aname: &str, pname: &str) -> Result<Vec<f64>, LagritError> {
    let buffer = mmfindbk(8, 1, aname, pname)?;
    let buffer_f64 =
        unsafe { std::slice::from_raw_parts(buffer.as_ptr() as *const f64, buffer.len() / 8) };
    Ok(buffer_f64.to_vec())
}

pub fn mmfindbk_string(aname: &str, pname: &str) -> Result<Vec<String>, LagritError> {
    let buffer = mmfindbk(1, STR_BUF_LEN, aname, pname)?;

    let buffer_str = unsafe {
        std::str::from_utf8_unchecked(std::slice::from_raw_parts(buffer.as_ptr(), buffer.len()))
    }
    .to_string();

    Ok(buffer_str
        .as_bytes()
        .chunks(STR_BUF_LEN)
        .map(|b| unsafe { std::str::from_utf8_unchecked(b) })
        .map(|s| s.trim().to_string())
        .collect::<Vec<String>>())
}

pub fn mmrelprt(pname: &str) -> Result<(), LagritError> {
    let pname = CString::new(pname)?;
    let mut status = 0;

    unsafe {
        fc_mmrelprt(pname.as_ptr(), &mut status);
    }

    if status != 0 {
        Err(status.into())
    } else {
        Ok(())
    }
}

pub fn cmo_get_index(mo: &str) -> Result<usize, LagritError> {
    let mo = CString::new(mo)?;
    let mut idx = -1;
    let mut status = 0;

    unsafe { fc_cmo_get_index(mo.as_ptr(), &mut idx, &mut status) };

    if status != 0 {
        return Err(status.into());
    } else if idx <= 0 {
        return Err(LagritError::MeshObjectNotFound);
    }

    Ok(idx as usize)
}

pub fn cmo_get_name() -> Result<String, LagritError> {
    let name = [0i8; STR_BUF_LEN];
    let mut status: i32 = 0;

    unsafe { fc_cmo_get_name(name.as_ptr() as *mut i8, &mut status) };

    if status != 0 {
        return Err(status.into());
    }

    let mo_name = unsafe {
        std::str::from_utf8_unchecked(std::mem::transmute::<&[i8], &[u8]>(name.as_slice()))
    };

    Ok(mo_name.trim_matches('\0').to_string())
}

pub fn cmo_get_mesh_type(mo: &str) -> Result<(String, i32), LagritError> {
    let mo = CString::new(mo)?;
    let mut mesh_type = [0i8; STR_BUF_LEN];
    let mut imesh_type = 0;
    let mut status = 0;

    unsafe {
        fc_cmo_get_mesh_type(
            mo.as_ptr(),
            mesh_type.as_mut_ptr(),
            &mut imesh_type,
            &mut status,
        )
    };

    if status != 0 {
        return Err(status.into());
    }

    let mesh_type = unsafe {
        std::str::from_utf8_unchecked(std::mem::transmute::<&[i8], &[u8]>(mesh_type.as_slice()))
    };

    Ok((mesh_type.trim_matches('\0').to_string(), imesh_type))
}

pub fn cmo_attlist(mo: &str) -> Result<Vec<AttrInfo>, LagritError> {
    Ok(mmfindbk_string("cmo_attlist", mo)?
        .chunks(7)
        .map(|c| AttrInfo {
            name: c[0].clone(),
            type_: c[1].clone(),
            rank: c[2].clone(),
            length: c[3].clone(),
            interpolation: c[4].clone(),
            persistence: c[5].clone(),
            ioflag: c[6].clone(),
        })
        .collect())
}

pub fn cmo_get_info(option: &str, mo: &str) -> Result<AttrValue, LagritError> {
    let mo = if mo == "-cmo-" || mo == "-default-" || mo == "-def-" {
        cmo_get_name()?
    } else {
        mo.to_string()
    };

    let _ = cmo_get_index(&mo)?;

    let option = match option {
        "imt" => "imt1",
        "icr" => "icr1",
        "itp" => "itp1",
        "isn" => "isn1",
        _ => option,
    };

    let mo_atts = cmo_attlist(&mo)?;

    for (idx, att) in mo_atts.iter().enumerate() {
        if att.name == option {
            if att.type_ == "INT" {
                return Ok(AttrValue::Int(
                    mmfindbk_i64("cmo_attparam_idefault", &mo)?[idx],
                ));
            } else if att.type_ == "VINT" {
                return Ok(AttrValue::Vint(mmfindbk_i64(option, &mo)?));
            } else if att.type_ == "REAL" {
                return Ok(AttrValue::Real(
                    mmfindbk_f64("cmo_attparam_rdefault", &mo)?[idx],
                ));
            } else if att.type_ == "VDOUBLE" {
                return Ok(AttrValue::Vdouble(mmfindbk_f64(option, &mo)?));
            } else if att.type_ == "CHARACTER" {
                return Ok(AttrValue::Character(
                    mmfindbk_string("cmo_attparam_cdefault", &mo)?[idx].clone(),
                ));
            } else if att.type_ == "VCHAR" {
                return Ok(AttrValue::Vchar(mmfindbk_string(option, &mo)?));
            }
        }
    }

    Err(LagritError::Unknown)
}

pub fn cmo_set_info(option: &str, mo: &str, value: AttrValue) -> Result<(), LagritError> {
    let attr = match option {
        "imt" => "imt1",
        "icr" => "icr1",
        "itp" => "itp1",
        "isn" => "isn1",
        _ => option,
    };

    let attr_info = cmo_attlist(mo)?
        .into_iter()
        .find(|info| info.name == attr)
        .ok_or(LagritError::AttributeNotFound)?;

    let attr_rank = match attr_info.rank.as_str() {
        "scalar" => 1,
        lenvar => match cmo_get_info(lenvar, mo)? {
            AttrValue::Int(v) => v,
            _ => return Err(LagritError::Unknown),
        },
    };

    let attr_len = match attr_info.length.as_str() {
        "scalar" => 1,
        lenvar => match cmo_get_info(lenvar, mo)? {
            AttrValue::Int(v) => v,
            _ => return Err(LagritError::Unknown),
        },
    };

    let attr = CString::new(attr)?;
    let mo = CString::new(mo)?;

    let mut status = 0;

    match value {
        AttrValue::Vint(mut v) => unsafe {
            if attr_info.type_ != "VINT" {
                return Err(LagritError::InvalidValueType);
            }
            let mut len = v.len() as i64;
            if len != attr_rank * attr_len {
                return Err(LagritError::InvalidValueLength);
            }
            fc_set_iattr(
                attr.as_ptr(),
                mo.as_ptr(),
                v.as_mut_ptr(),
                &mut len,
                &mut status,
            );
        },
        AttrValue::Vdouble(mut v) => unsafe {
            if attr_info.type_ != "VDOUBLE" {
                return Err(LagritError::InvalidValueType);
            }
            let mut len = v.len() as i64;
            if len != attr_rank * attr_len {
                return Err(LagritError::InvalidValueLength);
            }

            fc_set_fattr(
                attr.as_ptr(),
                mo.as_ptr(),
                v.as_mut_ptr(),
                &mut len,
                &mut status,
            );
        },
        _ => {
            return Err(LagritError::UnsupportedDataType);
        }
    }

    if status != 0 {
        return Err(status.into());
    }

    Ok(())
}
