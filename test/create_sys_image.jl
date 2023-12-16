@info "Loading packages ..."
using ModelingToolkit, OrdinaryDiffEq, PackageCompiler, PyPlot, Timers

@info "Creating sysimage ..."
push!(LOAD_PATH,joinpath(pwd(),"src"))

PackageCompiler.create_sysimage(
    [:ModelingToolkit, :OrdinaryDiffEq, :Timers];
    sysimage_path="kps-image_tmp.so",
    precompile_execution_file=joinpath("test", "test_for_precompile.jl")
)