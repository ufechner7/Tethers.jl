### Unreleased
#### Fixed
- the `run_julia` now modifies the LD_LIBRARY_PATH to avoid problems with shared libraries on Linux
- out-of memory error when running `create_sys_image` on systems with 16GB RAM
- error on Windows when using the `Tether_6c.jl` example
#### Changed
- renamed ODESystem to System
- updated the `install` script to also compile the docs project
- update the `create_sys_image` script

### Tethers v1.2.2 2026-03-20
#### Added
- add CITATION.cff

### Tethers v1.2.1 2026-03-20
#### Changed
- works now with Julia 1.11 or Julia 1.12; use `juliaup default 1.11` or `juliaup default 1.12` to select your preferred Julia version
- uses MTK 11, which is much, much faster when simplifying complex equation systems
- new `bin/install` script. Use it before running `bin/create_sys_image`.
- updated many other packages to the latest version
- the `install` and the `create_sys_image` now support the parameter `--yes` for non-interactive use

### Tethers v1.2.0 2025-10-22
#### Changed
- works now with Julia 1.10 or Julia 1.11; use `juliaup default 1.10` or `juliaup default 1.11` to select your preferred Julia version
- updated the package versions 