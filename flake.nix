{
  inputs = {
    fenix = {
      url = "github:nix-community/fenix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    flake-utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
    rust-overlay.url = "github:oxalica/rust-overlay";
  };

  outputs =
    { self
    , fenix
    , flake-utils
    , nixpkgs
    , rust-overlay
    ,
    }:
    flake-utils.lib.eachDefaultSystem (system:
    let
      toolchain = with fenix.packages.${system};
        combine [
          stable.cargo
          stable.clippy
          stable.rust-src
          stable.rustc
          stable.rustfmt
          rust-analyzer
        ];
      overlays = [ (import rust-overlay) ];
      pkgs = import nixpkgs {
        inherit system overlays;
      };
    in
    {
      devShells.default = pkgs.mkShellNoCC {
        buildInputs = [
          toolchain
          pkgs.cargo-audit
          pkgs.libclang.lib
          pkgs.llvmPackages.libcxxClang
          pkgs.cmake
          pkgs.gfortran
          pkgs.libcxx
          pkgs.libiconv
          pkgs.curl
          pkgs.maturin
          pkgs.python3Full
          pkgs.python3Packages.pip
          pkgs.python3Packages.numpy
          pkgs.python3Packages.pyvista
          pkgs.python3Packages.pytest
          pkgs.python3Packages.venvShellHook
        ];
        venvDir = "venv";
        env = {
          RUSTFLAGS = "-C linker=${pkgs.gfortran}/bin/gcc";
          CC = "${pkgs.gfortran}/bin/gcc";
          CXX = "${pkgs.gfortran}/bin/g++";
          LIBCLANG_PATH = "${pkgs.libclang.lib}/lib";
        };
      };
    });
}
