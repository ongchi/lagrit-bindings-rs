use std::env;

use lagrit_bindings::{InitMode, LaGriT, LagritError};

fn main() -> Result<(), LagritError> {
    let workdir = format!("{}/examples/ex1_output", env!["CARGO_MANIFEST_DIR"]);
    std::fs::create_dir_all(&workdir).unwrap();
    env::set_current_dir(workdir).unwrap();

    let lg = LaGriT::new(InitMode::Slient, None, None)?;

    lg.sendcmd("cmo/create/mocyl/ / /tet")?;
    lg.sendcmd("createpts/rtz/3,4,2/0.0,0.0,0.0/6.0,360.0,2.0/1,1,1")?;
    lg.sendcmd("pset/p1/attribute/xic/1,0,0/get 1.0")?;
    lg.sendcmd("pset/p2/attribute/xic/1,0,0/get 1.0")?;
    lg.sendcmd("pset/p3/attribute/xic/1,0,0/get 1.0")?;
    lg.sendcmd("connect/noadd")?;
    lg.sendcmd("dump/avs/test.inp/mocyl")?;

    let cmo = lg.cmo()?;
    println!("cmo: {}", cmo.name());
    cmo.status()?;

    println!("nnodes: {:?}", cmo.attr("nnodes")?); // INT
    println!("itettyp: {:?}", cmo.attr("itettyp")?); // VINT

    println!("xmax: {:?}", cmo.attr("xmax")?); // REAL
    println!("xic: {:?}", cmo.attr("xic")?); // VDOUBLE

    println!("geom_name: {:?}", cmo.attr("geom_name")?); // CHARACTER
    println!("psetnames: {:?}", cmo.attr("psetnames")?); // VCHAR

    Ok(())
}
