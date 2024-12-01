# LaGriT Bindings

Rust bindings for [LaGriT](https://lagrit.lanl.gov/) (Los Alamos Grid Toolbox).

This crate provides both drop-in replacements for PyLaGriT and a Rust API for LaGriT.

# Installation

## From Source

### Prerequisites

- cmake
- gfortran
- gcc
- rust
- python

```bash
git clone --recurse-submodules https://github.com/ongchi/lagrit-bindings-rs.git
export RUSTFLAGS="-C linker=gcc"
export CC=gcc
maturin build --release
pip install target/wheels/lagrit_bindings-*.whl
```
