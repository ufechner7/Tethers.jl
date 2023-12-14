# References

## Python
[Assimulo](https://jmodelica.org/assimulo/) offers 14 solvers
with good documentation for explicit and implicit problems.

## Julia
- [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) offers a unified interface to about 150 different solvers for a large range of problems. It wraps many existing open-source and commercial solvers, that have been implemented in C++ or Fortran and adds a growing number of native Julia solvers, many of them state-of-the-art.
- [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) is an acausal modeling framework for automatically parallelized scientific machine learning (SciML) in Julia. A computer algebra system for integrated symbolics for physics-informed machine learning and automated transformations of differential equations. 
- [KiteModels.jl](https://ufechner7.github.io/KiteModels.jl/stable/) implements kite models, connected to a tether for airborne wind energy applications. It uses the same algorithms as this tutorial, but it is not (yet) using ModelingToolkit. 
- [Working with Julia projects](https://ufechner7.github.io/2022/08/16/julia-projects.html) A must-read before creating your first project.

## Scientific papers
- Uwe Fechner, Rolf van der Vlugt, Edwin Schreuder, Roland Schmehl. (2015). [*Dynamic Model of a Pumping Kite Power System*](https://www.sciencedirect.com/science/article/pii/S0960148115003080)  describes the tether model used in this tutorial, but also a model of a complete kite power system with experimental validation. Renewable Energy. [Preprint](http://arxiv.org/abs/1406.6218).
- Yingbo Ma, Shashi Gowda, Ranjan Anantharaman, Chris Laughman, Viral Shah, and Chris Rackauckas. (2021). [*ModelingToolkit: A Composable Graph Transformation System For Equation-Based Modeling.*](https://arxiv.org/abs/2103.05244)
- Rackauckas, Christopher and Nie, Qing (2017). [*DifferentialEquations.jl--a performant and feature-rich ecosystem for solving differential equations in Julia}*](https://openresearchsoftware.metajnl.com/articles/10.5334/jors.151) Journal of Open Research Software.
- D.F. Duda1, H. Fuest, T. Islam, T. Ostermann, D. Moormann1. (2022). [*Hybrid modeling approach for the tether of an airborne wind energy system*](https://link.springer.com/article/10.1007/s13272-022-00581-7) CEAS Aeronautical Journal.