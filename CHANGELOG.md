### Tethers v2.1.0 (unreleased)

#### Added

- `examples/python/bench_casadi.py`, deriving the Jacobian of the `Tether_08`
  model with `ca.jacobian` instead of by hand, checking it against the
  hand-derived `calc_acc_jac` of `Tether_08.py`, and timing the SUNDIALS
  integrators CasADi ships; needs `pip install casadi` and is not part of the
  tutorial
- `examples/quasisteady/benchmark_scaling.jl` and
  `docs/images/qsm_vs_dynamic.png`, comparing the cost of the quasi-steady model
  against the dynamic one over 4 to 32 tether segments, both per simulated
  second; the quasi-steady solve is 5.6x cheaper at four segments and 49x at
  thirty-two
- `.github/workflows/draft-paper.yml`, building the JOSS paper with the
  journal's own `openjournals/inara` image on every push that touches `paper/`
  and uploading it as the `paper` artifact, so no local Docker or LaTeX
  toolchain is needed
- Bart van de Lint as a third author of the JOSS paper, and Andrea Bertozzi to
  `.zenodo.json`, which listed only one creator

#### Changed

- every example's time simulation now builds its `ODEProblem` with
  `jac=true, sparse=true`, so ModelingToolkit generates and compiles the
  analytic, sparse Jacobian ahead of time instead of the solver rebuilding a
  dense one by automatic differentiation at every step; the `autodiff=` keyword
  these solves passed to `FBDF` is now redundant and was dropped, and with it
  the `ADTypes` import where nothing else used it
- `examples/Tether_11.jl` keeps its finite-difference Jacobian: with a forward-mode or an analytic one, `FBDF` drives dt below eps at t = 0 and returns `Unstable`. This predates the change - plain `FBDF()` with no Jacobian at all already fails there - so the example is left as it was
- the steady-state solves are unchanged; only the time simulations were touched
- the Python examples no longer use Assimulo. Each model is now written once as
  a CasADi expression graph, CasADi derives its exact Jacobian and sparsity, and
  SUNDIALS' CVODES integrates it with a sparse Newton solve - the same
  fixed-leading-coefficient BDF family as the Julia examples' `FBDF`. The
  hand-written Jacobians of `Tether_06.py`, `Tether_06c.py`, `Tether_07.py` and
  `Tether_08.py`, 86 to 150 lines of calculus each, are gone; the ten examples
  lost about 700 lines between them. All ten still reproduce the Julia
  trajectories within the tolerances `test/test_tether_*.jl` apply
- `CondaPkg.toml` swaps `assimulo` for `casadi` and adds `scipy`, which
  `Tether_08.py`'s steady-state solve uses
- `examples/python/bench_casadi.py` now imports the model from `Tether_08.py`
  instead of defining its own copy
- `Tether_03b.py` locates the taut/slack crossings itself, by bisection between
  output points: CasADi's event detection (a `zero` entry in the DAE dictionary)
  is experimental and aborts on this model with "tout too far back in direction
  of integration". `Tether_06c.py`, whose events fire once per segment, uses it
  successfully
- `paper/paper.md` describes what the package now does: the Python examples use
  CasADi and SUNDIALS rather than Assimulo's IDA, both implementations use an
  analytic sparse Jacobian, and the performance comparison reports them within
  about 20% of each other instead of the earlier 13-30x, with an explanation of
  where that number came from. The AI usage disclosure covers the Jacobian and
  CasADi work
- `paper/paper.bib` cites CasADi [Andersson2019] instead of Assimulo, which the
  examples no longer use
- `paper/build` used the third-party `openbases/openbases-pdf` image, which does
  not produce the journal's layout; it now uses `openjournals/inara`, and points
  at the CI workflow as the primary route
- `docs/julia_vs_python.md` now compares like with like: both sides get an
  analytic sparse Jacobian and a BDF integrator, which closes most of the gap it
  used to report

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
