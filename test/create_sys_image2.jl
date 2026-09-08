# Build a system image from the active project. This variant is used when `Tethers` was
# installed with the package manager and the examples were installed with
# `install_examples()`; `copy_bin()` installs it as `test/create_sys_image.jl`.
#
# Unlike the script of the same name in a clone of the Tethers.jl repository, this script
# does not activate the `test` project; it uses the project that is already active, which
# is where `install_examples()` added the packages of the example scripts. `PackageCompiler`
# is not one of them, so it is added here if needed.
using Pkg
if ! ("PackageCompiler" ∈ keys(Pkg.project().dependencies))
    @info "Installing PackageCompiler ..."
    Pkg.add("PackageCompiler")
end

@info "Loading packages ..."
using ModelingToolkit, OrdinaryDiffEqCore, OrdinaryDiffEqBDF
using SteadyStateDiffEq, PackageCompiler, MakieControlPlots, Timers, REPL.TerminalMenus

FAST=true

@info "Creating sysimage ..."
push!(LOAD_PATH,joinpath(pwd(),"src"))

pkgs=[:ModelingToolkit, :OrdinaryDiffEqCore, :OrdinaryDiffEqBDF,
      :SteadyStateDiffEq, :Timers]
if FAST
    push!(pkgs, :MakieControlPlots)
end

GC.gc(true)
let mem = Sys.free_memory() / 1024^2
    @info "Free memory: $(round(mem; digits=1)) MB"
    if haskey(ENV, "JULIA_IMAGE_THREADS")
        @info "JULIA_IMAGE_THREADS: $(ENV["JULIA_IMAGE_THREADS"])"
    else
        @info "JULIA_IMAGE_THREADS not defined!"
    end
end

PackageCompiler.create_sysimage(
    pkgs;
    sysimage_path="kps-image_tmp.so",
    precompile_execution_file=joinpath("test", "test_for_precompile.jl")
)
