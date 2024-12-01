# LaGriT Bindings

Rust bindings for [LaGriT](https://lagrit.lanl.gov/) (Los Alamos Grid Toolbox).

This crate provides both drop-in replacements for PyLaGriT and a Rust API for LaGriT.

# Installation

```bash
git clone --recurse-submodules https://github.com/ongchi/lagrit-bindings-rs.git
maturin build --release
pip install target/wheels/lagrit_bindings-*.whl
```
