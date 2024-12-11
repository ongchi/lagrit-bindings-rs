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

### Clone from GitHub repository

```bash
git clone --recurse-submodules https://github.com/ongchi/lagrit-bindings-rs.git
```

```bash
cd lagrit-bindings-rs/lagrti-sys/vendor/lagrit
./install-exodus.sh
# For MacOS
# ./MAC_install-exodus.sh
```

### Build LaGriT Bindings

```bash
# (Optional) Change default linker to gcc on MacOS
export RUSTFLAGS="-C linker=gcc"
export CC=gcc
```

then build the bindings.

```bash
cd lagrit-bindings-rs
maturin build --release
pip install target/wheels/lagrit_bindings-*.whl
pip install -r requirements.in
```

### (Optional) Build with ExodusII support

```bash
cd lagrit-bindings-rs/lagrit-sys/vendor/lagrit
./install-exodus.sh
# For MasOS
# ./MAC_install-exodus.sh
```

then build the bindings with exodus feature.

```bash
maturin build --release --features exodus
```

### Install LaGriT Binary Executable

```bash
cargo install --example lagrit --features=clap,exodus --path .
```
