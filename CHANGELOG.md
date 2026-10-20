### Tethers unreleased
#### Added
- `bin/mem_monitor`, logging RAM, swap, zram and CPU usage while a command runs and printing its peak memory demand; it replaces `bin/peak_mem`
- `bin/create_sys_image_logged`, building the system image under `mem_monitor` and saving the output, the peak memory demand and the total run time in `sysimage.log`
- `.JETLSConfig.toml.default`, the JETLS configuration; `bin/install` copies it to `.JETLSConfig.toml` if that does not exist

#### Changed
- the number of threads for compiling the system image no longer grows with the memory: measurements showed that more than two threads need no extra memory, so with at least 35 GB of RAM + swap (including the size the pagefile can grow to on Windows) half of the CPU threads are used, at most 16, otherwise one; on Linux and macOS, machines with 16 to 27 GB of RAM use two threads instead of one (one if they have fewer than 4 CPU threads), and only machines with less than 16 GB of RAM are limited to one thread
- `create_sys_image` stops with an error on all platforms, including Windows, if there is less than 30 GB of RAM + swap, which is not enough even for one thread; on Linux the bash script checks this before doing anything else
- `create_sys_image` respects `JULIA_IMAGE_THREADS` when it is set by the caller
- `create_sys_image` runs all `julia` calls with `--startup-file=no`
- CI prints the CPU and compiles for the portable CPU target `x86-64-v3` on x64, so a package cache built on a newer runner CPU is no longer rejected on an older one
- fixed JETLS diagnostics in `src` and the examples

#### Fixed
- `create_sys_image` failed on macOS while checking the memory, because it read `/proc/meminfo`; there it now counts the RAM and the free space of the volume that macOS creates its swap files on

### Tethers v2.1.0 2026-10-18
#### Added
- `bin/install` options `--update` (regenerate the manifest with `Pkg.update()` instead of instantiating the pinned one), `--test` (run the tests after installing) and `--version <ver>`/`+<ver>` (use that Julia version without changing the juliaup default)
- `bin/update_default_manifests`, regenerating `Manifest-v<ver>.toml.default` for each supported Julia version, only if the install and the tests succeed
- documentation of `Tether_11.jl` in the examples section of the docs

#### Changed
- dropped support for Julia 1.11; Julia 1.12 and 1.13 are supported
- the system image is built for the native CPU on Julia 1.13 too; the PackageCompiler default there compiles every function for several CPU targets, needing more than twice the memory and build time
- building the system image works now on Windows with Julia 1.13
- CI tests Julia 1.12 and 1.13 on Linux, and Julia 1.13 on Windows and macOS
- regenerated the pinned `Manifest-v1.12.toml.default` and `Manifest-v1.13.toml.default`

### Tethers v2.0.0 2026-09-11
#### Added
- `examples/Project.toml` and `test/Project.toml`, joined to the root package as workspace members on Julia 1.12, so the example and test dependencies no longer bloat the main `Project.toml`
- `example_packages()` in `Tethers.jl`, deriving the package list for `install_examples` from `examples/Project.toml` instead of a hardcoded (and, until now, incomplete) list
- `Tethers.QuasiSteady` submodule (`src/Tether_quasisteady.jl`), a quasi-steady tether model based on Williams (2017); exported API `StaticSettings`, `Tether`, `init!`, `step!`
- `src/qsm_conventions.jl`, converting between the MATLAB reference fixtures' angle convention and this package's [KiteUtils.jl reference frames](https://opensourceawe.github.io/KiteUtils.jl/stable/reference_frames/) convention
- `examples/quasisteady/` (`benchmark_qsm.jl`, `flying_circular.jl`, `force_plots.jl`, `run_catenary.jl`, `run_catenary_matlab.jl`) and `examples/menu3.jl` to run them
- `examples/Tether_11.jl`, a kite flying a circular trajectory on a cone with the dynamic mass-spring-damper tether, for comparison against `examples/quasisteady/flying_circular.jl`; added to `menu.jl`
- `MovingEnd` and `assemble_tether` in `TetherComponent.jl`, imposing a prescribed, time-dependent trajectory on a tether end
- `test/test_qsm.jl` and `test/data/*.mat`, 70 assertions checking the quasi-steady model against MATLAB reference data
- `docs/quasisteady.md`, documenting the port of the quasi-steady model, the investigations along the way and their outcomes
- Julia 1.13 support, with the pinned `Manifest-v1.13.toml.default`
- reference to Williams, *Cable Modeling Approximations for Rapid Simulation* (2017), in `README.md` and `docs/src/references.md`

#### Fixed
- `install_examples` no longer omits `Symbolics`, which several examples need directly
- `copy_examples` no longer copies `examples/Project.toml` or leftover `Manifest*.toml` files into the destination directory
- angle convention mismatch between `res!` and the MATLAB reference fixtures (a parametrisation difference, not a frame rotation); the four reference comparisons in `test/test_qsm.jl`, previously `@test_broken`, now pass at `rtol=1e-9`
- the `.mat` fixtures' `T.rho_t` is a mass per unit length, not a density as `Settings.rho_tether` is, so the tether had been simulated about 1442 times too light
- `maxiters` warning in `Tether_11.jl`, caused by two independent bugs; `Tether_09.jl` and `Tether_11.jl` now check the steady-state solve's return code instead of silently reusing a non-converged state
- `Tether_11.jl` was dead code (its body sat inside a triple-quoted string) and did not build under MTK 11 (`vcat` of heterogeneous equations produced an uninferable `Vector{Any}`)
- stale `PreallocationTools` compat bound that made the workspace root unresolvable against 1.x
- CI no longer uses `julia-actions/julia-runtest`, whose `Pkg.test()` sandbox broke the `[sources]` path in `test/Project.toml`; the test project is now activated and run in place

#### Breaking
- the main `Project.toml` now only lists the packages needed by `src/`; `ADTypes`, `GLMakie`, `LaTeXStrings`, `LinearSolve`, `MakieControlPlots`, `OrdinaryDiffEq`, `PackageCompiler`, `StatsBase`, `SteadyStateDiffEq`, `Test` and `Timers` moved to `examples/Project.toml` and/or `test/Project.toml`

#### Changed
- `LiveServer` is no longer a dependency of Tethers; `docu()` now expects it in your default Julia environment, and `bin/install` adds it there if missing
- `bin/run_julia` and `bin/create_sys_image` now activate `examples`/`test` instead of the root project, matching where their dependencies now live
- `bin/install` also instantiates the `examples` and `test` subprojects; on Julia 1.11 (whose `Pkg` does not support workspaces) `docs/` is seeded from the pinned root manifest as before, while `examples/` and `test/` resolve independently since seeding them from the now much smaller root manifest can make their resolve unsatisfiable
- regenerated the pinned `Manifest-v1.11.toml.default` and `Manifest-v1.12.toml.default` to match the new project structure
- `simulate_tether` is about 4x faster (23.4 µs vs 94 µs, 63 vs 1118 allocations, 22 vs 36 solver iterations), from a non-allocating residual, solving for tension on a logarithmic scale with a capped trust region, and deriving the initial tension guess from the catenary fit instead of a fixed fraction of the spring constant
- removed `src/Tether_qsm_dual.jl`, the dual-number copy of the quasi-steady model, now that the main model differentiates cleanly
- Julia 1.13 is now the default; CI caches the whole Julia depot instead of just artifacts, and excludes it from Windows Defender scanning to speed up precompilation on Windows

### Tethers v1.2.3 2026-09-08
#### Added
- `TetherComponent.jl`, a re-usable tether component (`Tether`, `FixedEnd`, `FreeEnd`, `TetherSettings`, `set_diameter!`, `m_end`, `Point3D`), exported from `Tethers.TetherComponents` and used by the new `Tether_10.jl` example
- new Python examples `Tether_02.py`, `Tether_06c.py` and `Tether_08.py`, matching their Julia counterparts
- hand-coded, analytic Jacobians for the Python examples, replacing finite-difference Jacobians
- test scripts for all tether examples (`test_tether_01.jl` .. `test_tether_08.jl`, `test_tether_10.jl`) that check that the Python and Julia implementations produce the same results
- `test/test_copy_install.jl`, unit tests for `copy_files`, `copy_examples`, `copy_bin` and `install_examples`
- `test/test_tether_component.jl`, unit tests for the re-usable tether component, checking its steady state, drag, catenary shape and compression stiffness against analytic results
- `docs/julia_vs_python.md`, comparing the performance and code size of the Julia and Python implementations
- `install_examples`, `copy_examples` and `copy_bin` functions in `Tethers.jl`, to install the example scripts and helper scripts (`bin/run_julia`, `bin/install`, `bin/create_sys_image`) into the current working directory, and optionally add the packages they need
- `examples/menu.jl` and `examples/menu2.jl`, replacing `src/init.jl` as the entry points for running the examples interactively
- code coverage reporting to CI, uploaded to Codecov
#### Fixed
- out-of memory error when running `create_sys_image` on systems with 16GB RAM
- error on Windows when using the `Tether_6c.jl` example
- numerous bugs in the Python examples (`Tether_01.py` .. `Tether_08.py`), now producing results consistent with the Julia versions

#### Changed
- renamed ODESystem to System
- update the `create_sys_image` script; the GC heap size hint now scales with the available RAM instead of always being 8000M
- switched the interactive plots from PyPlot/matplotlib to MakieControlPlots; removed the PyCall and Conda dependencies from `Project.toml` and the Conda/matplotlib setup from `bin/run_julia`
- greatly improved the performance of the Python examples by using analytic Jacobians instead of finite-difference ones; Python is now only 13-30 times slower than Julia (15-20x on average), down from a much larger gap
- moved the Julia example scripts (`Tether_01.jl` .. `Tether_10.jl`) and the Python example scripts (`Tether_01.py` .. `Tether_08.py`, now under `examples/python`) from `src` to `examples`
- removed `src/init.jl` and the boilerplate at the top of the `RunTether_*.jl` scripts

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
