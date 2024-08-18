@info "Loading packages ..."
using ModelingToolkit, OrdinaryDiffEq, SteadyStateDiffEq, PackageCompiler, ControlPlots, Timers, REPL.TerminalMenus

FAST=true

@info "Creating sysimage ..."
push!(LOAD_PATH,joinpath(pwd(),"src"))

pkgs=[:ModelingToolkit, :OrdinaryDiffEq, :SteadyStateDiffEq, :Timers]
if FAST
    push!(pkgs, :ControlPlots)
end 

PackageCompiler.create_sysimage(
    pkgs;
    sysimage_path="kps-image_tmp.so",
    precompile_execution_file=joinpath("test", "test_for_precompile.jl")
)