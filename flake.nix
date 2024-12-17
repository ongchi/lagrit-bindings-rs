{
  inputs = {
    fenix = {
      url = "github:nix-community/fenix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    flake-utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-24.11-darwin";
  };

  outputs = {
    self,
    fenix,
    flake-utils,
    nixpkgs,
  }:
    flake-utils.lib.eachDefaultSystem (system: let
      toolchain = fenix.packages.${system}.default.toolchain;
      pkgs = nixpkgs.legacyPackages.${system};
    in {
      devShells.default = pkgs.mkShellNoCC {
        buildInputs = [
          toolchain
          pkgs.cmake
          pkgs.gfortran
          pkgs.libcxx
          pkgs.libiconv
          pkgs.curl
          pkgs.maturin
          pkgs.python3Packages.venvShellHook
        ];
        venvDir = "venv";
      };
    });
}
