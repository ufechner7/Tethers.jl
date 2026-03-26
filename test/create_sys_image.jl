@info "Loading packages ..."
using ModelingToolkit, OrdinaryDiffEq, SteadyStateDiffEq, PackageCompiler, ControlPlots, Timers, REPL.TerminalMenus

FAST=true

@info "Creating sysimage ..."
push!(LOAD_PATH,joinpath(pwd(),"src"))

pkgs=[:ModelingToolkit, :OrdinaryDiffEq, :SteadyStateDiffEq, :Timers]
if FAST
    push!(pkgs, :ControlPlots)
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