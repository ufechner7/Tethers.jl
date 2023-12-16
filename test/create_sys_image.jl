@info "Loading packages ..."
using ModelingToolkit, OrdinaryDiffEq, PackageCompiler, PyPlot, Timers

FAST=true

@info "Creating sysimage ..."
push!(LOAD_PATH,joinpath(pwd(),"src"))

pkgs=[:ModelingToolkit, :OrdinaryDiffEq, :Timers]
if FAST
    push!(pkgs, :PyPlot)
end 

PackageCompiler.create_sysimage(
    pkgs;
    sysimage_path="kps-image_tmp.so",
    precompile_execution_file=joinpath("test", "test_for_precompile.jl")
)