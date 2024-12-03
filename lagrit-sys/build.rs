use std::env;
use std::path::PathBuf;

fn main() {
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    let vendor_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap()).join("vendor");

    // LaGriT
    let mut cmake_config = cmake::Config::new(vendor_dir.join("lagrit"));
    cmake_config.define("CMAKE_INSTALL_PREFIX", &out_path);

    if cfg!(feature = "exodus") {
        cmake_config.define("LAGRIT_BUILD_EXODUS", "ON");
    }

    cmake_config.out_dir(out_path.join("build"));
    let lagrit_root = cmake_config.build();

    // Fortran bindings
    println!("cargo:rerun-if-changed=src/bindings.h");
    println!("cargo:rerun-if-changed=src/bindings.f90");

    cc::Build::new()
        .file("src/bindings.f90")
        .include("vendor/lagrit/src")
        .flag("-fcray-pointer")
        .compile("lagrit-bindings");

    let bindings = bindgen::Builder::default()
        .clang_arg(format!("-I{}", "vendor/lagrit/src"))
        .header("src/bindings.h")
        .generate()
        .expect("Unable to generate bindings.");

    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Unable to write bindings.");

    // Linking
    println!(
        "cargo:rustc-link-search=native={}",
        lagrit_root.join("build").display()
    );
    println!("cargo:rustc-link-lib=static=lagrit");
    println!("cargo:rustc-link-lib=static=lagrit-bindings");
    println!("cargo:rustc-link-lib=gfortran");
    println!("cargo:rustc-link-lib=c");
    println!("cargo:rustc-link-lib=c++");
    println!("cargo:rustc-link-lib=stdc++");

    if cfg!(feature = "exodus") {
        println!(
            "cargo:rustc-link-search=native={}",
            vendor_dir
                .join("lagrit")
                .join("TPLs")
                .join("seacas")
                .join("lib")
                .display()
        );
        println!("cargo:rustc-link-lib=static=exodus");
        println!("cargo:rustc-link-lib=static=exoIIv2for32");
        println!("cargo:rustc-link-lib=static=netcdf");
        println!("cargo:rustc-link-lib=static=hdf5");
        println!("cargo:rustc-link-lib=static=hdf5_hl");
    }
}
