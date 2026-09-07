### Tethers v1.2.3 2026-09-08
#### Added
- `TetherComponent.jl`, a re-usable tether component (`Tether`, `FixedEnd`, `FreeEnd`, `TetherSettings`, `set_diameter!`, `m_end`, `Point3D`), exported from `Tethers.TetherComponents` and used by the new `Tether_10.jl` example
- new Python examples `Tether_02.py`, `Tether_06c.py` and `Tether_08.py`, matching their Julia counterparts
- hand-coded, analytic Jacobians for the Python examples, replacing finite-difference Jacobians
- test scripts for all tether examples (`test_tether_01.jl` .. `test_tether_08.jl`, `test_tether_10.jl`) that check that the Python and Julia implementations produce the same results
- `test/test_copy_install.jl`, unit tests for `copy_files`, `copy_examples`, `copy_bin` and `install_examples`
- `bin/create_pdf`, which renders a markdown document to PDF with pandoc and xelatex
- `test/test_tether_component.jl`, unit tests for the re-usable tether component, checking its steady state, drag, catenary shape and compression stiffness against analytic results
- `docs/julia_vs_python.md`, comparing the performance and code size of the Julia and Python implementations
- `install_examples`, `copy_examples` and `copy_bin` functions in `Tethers.jl`, to install the example scripts and helper scripts (`bin/run_julia`, `bin/install`, `bin/create_sys_image`) into the current working directory, and optionally add the packages they need
- `examples/menu.jl` and `examples/menu2.jl`, replacing `src/init.jl` as the entry points for running the examples interactively
- code coverage reporting to CI, uploaded to Codecov
- `examples/Project.toml` and `test/Project.toml`, joined to the root package as workspace members on Julia 1.12, so the example and test dependencies no longer bloat the main `Project.toml`
- `example_packages()` in `Tethers.jl`, deriving the package list for `install_examples` from `examples/Project.toml` instead of a hardcoded (and, until now, incomplete) list
#### Fixed
- out-of memory error when running `create_sys_image` on systems with 16GB RAM
- error on Windows when using the `Tether_6c.jl` example
- numerous bugs in the Python examples (`Tether_01.py` .. `Tether_08.py`), now producing results consistent with the Julia versions
- `install_examples` no longer omits `Symbolics`, which several examples need directly
- `copy_examples` no longer copies `examples/Project.toml` or leftover `Manifest*.toml` files into the destination directory
#### Changed
- the position of `FixedEnd` is now the parameter `pos_fix` instead of a literal, so that a compiled model can be re-solved for a different anchor position without calling `mtkcompile` again
- the wind of `Tether` is now the parameter `v_wind` instead of a literal, so that a compiled model can be re-solved for a different wind speed without calling `mtkcompile` again
- the tether cross section of `Tether` is now given by the parameters `d_tether`, `c_spring_unit`, `damping_unit` and `mass_per_m` instead of literals, so that a compiled model can be re-solved for a different tether diameter without calling `mtkcompile` again
- renamed ODESystem to System
- update the `create_sys_image` script; the GC heap size hint now scales with the available RAM instead of always being 8000M
- switched the interactive plots from PyPlot/matplotlib to MakieControlPlots; removed the PyCall and Conda dependencies from `Project.toml` and the Conda/matplotlib setup from `bin/run_julia`
- greatly improved the performance of the Python examples by using analytic Jacobians instead of finite-difference ones; Python is now only 13-30 times slower than Julia (15-20x on average), down from a much larger gap
- moved the Julia example scripts (`Tether_01.jl` .. `Tether_10.jl`) and the Python example scripts (`Tether_01.py` .. `Tether_08.py`, now under `examples/python`) from `src` to `examples`
- removed `src/init.jl` and the boilerplate at the top of the `RunTether_*.jl` scripts
- the main `Project.toml` now only lists the packages needed by `src/`; `ADTypes`, `GLMakie`, `LaTeXStrings`, `LinearSolve`, `MakieControlPlots`, `OrdinaryDiffEq`, `PackageCompiler`, `StatsBase`, `SteadyStateDiffEq`, `Test` and `Timers` moved to `examples/Project.toml` and/or `test/Project.toml`
- `LiveServer` is no longer a dependency of Tethers; `docu()` now expects it in your default Julia environment, and `bin/install` adds it there if missing
- `bin/run_julia` and `bin/create_sys_image` now activate `examples`/`test` instead of the root project, matching where their dependencies now live
- `bin/install` also instantiates the `examples` and `test` subprojects; on Julia 1.11 (whose `Pkg` does not support workspaces) `docs/` is seeded from the pinned root manifest as before, while `examples/` and `test/` resolve independently since seeding them from the now much smaller root manifest can make their resolve unsatisfiable
- regenerated the pinned `Manifest-v1.11.toml.default` and `Manifest-v1.12.toml.default` to match the new project structure

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
