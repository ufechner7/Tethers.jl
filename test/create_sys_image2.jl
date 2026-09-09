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
    # Windows refuses to load a PE image of 2 GiB or more ("%1 is not a valid Win32
    # application"), and MakieControlPlots drags in both Makie backends. Listing GLMakie
    # instead keeps the interactive backend but leaves CairoMakie out of the image.
    push!(pkgs, Sys.iswindows() ? :GLMakie : :MakieControlPlots)
end

# Dropping docstrings and source-location metadata shrinks the image by a double-digit
# percentage, which Windows needs to stay under the 2 GiB limit.
build_args = Sys.iswindows() ? `--strip-metadata` : ``

function total_ram_swap_gb()
    if Sys.iswindows()
        # MEMORYSTATUSEX layout (x64): dwLength(u32), dwMemoryLoad(u32), then 7 x UInt64
        buf = zeros(UInt8, 64)
        unsafe_store!(Ptr{UInt32}(pointer(buf)), UInt32(64))
        ok = ccall((:GlobalMemoryStatusEx, "kernel32"), stdcall, Cint, (Ptr{UInt8},), buf)
        ok == 0 && error("GlobalMemoryStatusEx failed")
        total_pagefile = unsafe_load(Ptr{UInt64}(pointer(buf) + 24)) # ullTotalPageFile: RAM + pagefile
        total_pagefile / 1_073_741_824  # -> GiB
    else
        info = read("/proc/meminfo", String)
        memtotal = parse(Int, match(r"MemTotal:\s+(\d+)", info).captures[1])   # kB
        swaptotal = parse(Int, match(r"SwapTotal:\s+(\d+)", info).captures[1]) # kB
        (memtotal + swaptotal) / 1_048_576  # -> GiB
    end
end

let total = total_ram_swap_gb()
    @info "Total RAM + swap: $(round(total; digits=1)) GB"
    if total < 30
        msg = "At least 30 GB of RAM + swap is recommended to create a system image, " *
              "but only $(round(total; digits=1)) GB is available."
        if Sys.iswindows()
            # Windows' pagefile is commonly "System managed" and can grow on demand,
            # so a low reading here isn't necessarily a hard limit.
            @warn msg * " Windows may grow the pagefile automatically; increase it manually if the build fails."
        else
            error(msg * " Increase your swap file and retry.")
        end
    end

    if haskey(ENV, "JULIA_IMAGE_THREADS")
        @info "JULIA_IMAGE_THREADS already set to $(ENV["JULIA_IMAGE_THREADS"]), leaving as is"
    else
        # Linear interpolation: 1 thread at 30 GB, 8 threads at 40 GB (empirically the
        # amount of memory each extra image-compilation thread needs), capped at 16 and
        # at Julia's own default cap of CPU_THREADS / 2.
        ideal = 1 + floor(Int, (total - 30) * 7 / 10)
        threads = clamp(ideal, 1, 16)
        threads = min(threads, max(1, Sys.CPU_THREADS ÷ 2))
        ENV["JULIA_IMAGE_THREADS"] = string(threads)
        @info "Setting JULIA_IMAGE_THREADS=$threads based on $(round(total; digits=1)) GB total memory"
    end
end

GC.gc(true)
let mem = Sys.free_memory() / 1024^2
    @info "Free memory: $(round(mem; digits=1)) MB"
    @info "JULIA_IMAGE_THREADS: $(ENV["JULIA_IMAGE_THREADS"])"
end

PackageCompiler.create_sysimage(
    pkgs;
    sysimage_path="kps-image_tmp.so",
    precompile_execution_file=joinpath("test", "test_for_precompile.jl"),
    sysimage_build_args=build_args
)

let size_gib = filesize("kps-image_tmp.so") / 1024^3
    @info "System image size: $(round(size_gib; digits=2)) GiB"
    if Sys.iswindows() && size_gib > 1.8
        error("The system image is $(round(size_gib; digits=2)) GiB. Windows cannot load a " *
              "PE image of 2 GiB or more; it fails with \"%1 is not a valid Win32 application\". " *
              "Remove packages from `pkgs` or shorten test/test_for_precompile.jl.")
    end
end
