### Tethers v1.2.3 2026-09-06
#### Added
- `TetherComponent.jl`, a re-usable tether component (`Tether`, `FixedEnd`, `FreeEnd`, `TetherSettings`, `set_diameter!`, `m_end`, `Point3D`), exported from `Tethers.TetherComponents` and used by the new `Tether_10.jl` example
- new Python examples `Tether_02.py`, `Tether_06c.py` and `Tether_08.py`, matching their Julia counterparts
- hand-coded, analytic Jacobians for the Python examples, replacing finite-difference Jacobians
- test scripts for all tether examples (`test_tether_01.jl` .. `test_tether_08.jl`, `test_tether_10.jl`) that check that the Python and Julia implementations produce the same results
- `docs/julia_vs_python.md`, comparing the performance and code size of the Julia and Python implementations
#### Fixed
- out-of memory error when running `create_sys_image` on systems with 16GB RAM
- error on Windows when using the `Tether_6c.jl` example
- numerous bugs in the Python examples (`Tether_01.py` .. `Tether_08.py`), now producing results consistent with the Julia versions
#### Changed
- renamed ODESystem to System
- update the `create_sys_image` script
- switched the interactive plots from PyPlot/matplotlib to MakieControlPlots; removed the PyCall and Conda dependencies from `Project.toml` and the Conda/matplotlib setup from `bin/run_julia`
- greatly improved the performance of the Python examples by using analytic Jacobians instead of finite-difference ones; Python is now only 13-30 times slower than Julia (15-20x on average), down from a much larger gap

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