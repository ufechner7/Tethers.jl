---
title: 'Tethers.jl: A tutorial and a re-usable component for the simulation of tethers and cables'
tags:
  - Julia
  - tether
  - cable
  - mooring line
  - multibody dynamics
  - differential-algebraic equations
  - airborne wind energy
authors:
  - name: Uwe Fechner
    orcid: 0009-0008-2532-9458
    affiliation: 1
  - name: Andrea Bertozzi
    orcid: 0009-0005-1214-1737
    affiliation: 1, 2
  - name: Bart van de Lint
    affiliation: 3
affiliations:
  - name: Delft University of Technology, The Netherlands
    index: 1
  - name: Politecnico di Milano, Italy
    index: 2
  - name: Open Source AWE
    index: 3
date: 8 September 2026
bibliography: paper.bib
---

# Summary

Tethers, cables and ropes appear in many engineering problems: cranes, undersea cables,
the mooring lines of floating wind turbines, the launching of sailplanes and airborne wind
energy systems. Simulating them is harder than the simple picture of "a chain of masses
and springs" suggests. A modern high-modulus rope such as Dyneema is extremely stiff in
tension and almost incapable of carrying compression, which turns the equations of motion
into a stiff, non-smooth system that many general-purpose simulation environments handle
badly or not at all.

`Tethers.jl` addresses this in three ways. First, it is a tutorial: a graded series of
runnable examples that starts with a point mass thrown upwards and ends with a segmented,
reeling tether with aerodynamic drag and free end points. Each example is a short,
self-contained script, and most of them come with a matching Python implementation, so
that the modelling choices, the solver behaviour, the code size and the runtime can be
compared directly between the two ecosystems. Second, the package exports a re-usable,
acausal tether component built with ModelingToolkit.jl [@Ma2021], which can be connected to
other components of a larger model instead of being re-derived by hand for every new
application. Third, it provides a quasi-steady tether model after @Williams2017, which
replaces the time integration of the tether dynamics by a small nonlinear solve for the
tether shape and the forces at both ends, given the position and velocity of the far end.
It runs in tens of microseconds per step and is validated against the analytic catenary
and against reference results from the original MATLAB implementation. The dynamic and the
quasi-steady models share the same physical parameters, and an example that flies a kite on
a circular trajectory with both models allows their results to be compared directly.

![A five-segment tether, hanging from a fixed point and swinging under gravity, as
produced by one of the tutorial examples.\label{fig:tether}](tether_5segments.png){ width=70% }

# Statement of need

Tether models are a recurring building block in research on airborne wind energy, marine
systems and cable-driven machines, but they are usually written from scratch inside a
larger, application-specific code base and are rarely published on their own
[@Fechner2015; @Duda2022]. Researchers who need a tether therefore face two problems: there
is little material that explains *how* such a model is derived and why the naive
formulations fail, and there is no small, well-tested component they can simply attach to
their own model.

`Tethers.jl` was written to fill both gaps. The tutorial part makes the derivation explicit
and incremental. The model is a chain of $n$ spring-damper segments connecting $n+1$ point
masses. Reeling in or out changes the unstretched segment length, and with it the segment
stiffness, the segment damping and the particle mass, which is a common source of errors in
hand-written models. Slackness is represented by giving a compressed segment a much smaller
stiffness than a taut one, controlled by a relative compression stiffness parameter. Only
the component of the apparent wind perpendicular to a segment produces drag, and half of a
segment's drag is applied to each of its two particles. All of this is documented as
equations alongside the code, so the tutorial can be followed by readers who intend to
implement the model in a different language.

The component part packages the same physics for re-use. The connector `Point3D` carries
the position as its across variable and the force as its flow variable; a `Tether`
component with two such connectors is joined to a `FixedEnd` or a `FreeEnd`, and two
tethers can be chained by connecting both to the same `FreeEnd`. The wind velocity, the
anchor position and the tether cross section are model parameters rather than literals, so
a compiled model can be re-solved for a different wind speed, anchor position or tether
diameter without repeating the symbolic compilation step. The same algorithms, in a
hand-written form, are used by `KiteModels.jl` [@KiteModels] for airborne wind energy
simulations.

# Julia and Python side by side

Every tutorial example exists twice: as a Julia script using ModelingToolkit.jl and the
solvers of DifferentialEquations.jl [@Rackauckas2017; @Bezanson2017], and as a Python
script using CasADi [@Andersson2019] and the CVODES solver of SUNDIALS [@Hindmarsh2005]
that it ships with. Both sides integrate with a backward differentiation formula and both
derive an exact, sparse Jacobian from the model itself rather than by hand, and the test
suite checks that the two produce the same trajectories. This makes the comparison a
controlled one rather than an anecdote.

A stiff tether needs that Jacobian, and neither ecosystem supplies one by default. In
Julia, `ODEProblem(sys, ...; jac=true, sparse=true)` makes ModelingToolkit generate and
compile it ahead of time; without those two keywords the solver rebuilds a dense Jacobian
by automatic differentiation at every step, which costs a factor of 1.4 at five segments
and 7.2 at forty. In Python, `ca.jacobian` derives it from the model expression and finds
its sparsity. The Jacobian of the forty-segment model is 4% dense, and exploiting that
matters more than the analytic derivative itself: at forty segments the analytic Jacobian
alone gains about 15% over automatic differentiation, because a dense factorisation then
dominates the cost, while adding sparsity gains a factor of 7.2.

For a ten-second simulation of the ten-segment model, sampled every 20 ms at a relative and
absolute tolerance of $10^{-6}$, the two implementations take 17.1 ms and 20.7 ms; at forty
segments, 603 ms and 617 ms. They are within about 20% of each other in either direction,
and the differences between the ecosystems lie elsewhere. The Julia examples are roughly
half the length in lines of code and never write a Jacobian at all. Julia pays a one-time
compilation cost of seconds to tens of seconds, of which 0.4 to 4.1 seconds is the symbolic
Jacobian, so a script that solves once pays more than it saves, and installing the Julia
stack takes considerably longer than installing the Python one. Event handling favours
Julia: a `ContinuousCallback` is four lines, while CasADi's event support is experimental
and fails on one of the two examples that need it, which locates the taut/slack crossing by
bisection instead. Presenting both sides lets readers make an informed choice instead of
taking the benchmark on faith.

# Functionality

The package provides fifteen Julia examples and ten matching Python versions, reachable from an
interactive menu, plus five examples for the quasi-steady model in a second menu; the
`TetherComponents` submodule with the `Point3D`, `Tether`, `FixedEnd`, `FreeEnd` and
`MovingEnd` components and the `TetherSettings` parameter struct; the `QuasiSteady`
submodule with the `StaticSettings` parameter struct, the `Tether` state and the `init!`
and `step!` functions; helper functions to copy the examples and launcher scripts into a
user's own project; and a test suite that compares the Julia and Python results, checks
the component against analytic results for the steady state, the drag and the catenary
shape, and checks the quasi-steady model against the analytic catenary and MATLAB
reference data; and benchmarks for the Jacobian strategies of both ecosystems and for the
cost of the quasi-steady model against the dynamic one as the segment count grows.
Documentation, including the full derivation and all examples, is published online.

# AI usage disclosure

Generative AI tools were used in the development of this software and in the preparation
of this paper, as follows.

*Software.* The tutorial examples in their original form, the first Python implementations
and the re-usable tether component were written by the authors without AI assistance. The
quasi-steady model was first ported from MATLAB to Julia by hand. Claude Code (Anthropic)
was subsequently used as a coding assistant for four pieces of work: refactoring the
quasi-steady port to the exported `init!`/`step!` API and the performance work on its
nonlinear solve; giving the Julia examples an analytic sparse Jacobian and measuring what
that is worth; rewriting the Python examples around CasADi, which replaced the Jacobians
that had previously been derived by hand; and the benchmarks behind the comparison in the
section above. Every AI-assisted change was reviewed by the authors. The results were
verified by the test suite, which compares every Python example against its Julia
counterpart and the quasi-steady model against the analytic catenary and the MATLAB
reference results. GitHub Copilot's automated pull-request review was used to flag issues
in some changes; its suggestions were evaluated and applied by the authors.

*Documentation.* Parts of the documentation of the quasi-steady model, and the comparison
of the two ecosystems in `docs/julia_vs_python.md`, were drafted with Claude Code and
edited by the authors. The tutorial text and the derivation of the model were written by
the authors.

*Paper.* The authors wrote this paper. Claude Code was used to update the summary and
functionality sections after the quasi-steady model was added, to rewrite the comparison
section once both implementations used an analytic sparse Jacobian, to check citation
metadata for consistency, and to draft this disclosure. All text was reviewed and edited by
the authors, who take full responsibility for its content.

# Acknowledgements

The tether model implemented here originates in earlier work on kite power systems at Delft
University of Technology. The authors thank the developers of ModelingToolkit.jl and the
wider SciML ecosystem, on which this package depends.

This work has been supported by the MERIDIONAL project, which receives funding from the European Union’s Horizon Europe Program under the grant agreement no. [101084216](https://doi.org/10.3030/101084216). The opinions expressed in this document reflect only the author’s view and reflects in no way the European Commission’s opinions. The European Commission is not responsible for any use that may be made of the information it contains.

# References
