using Test
using Tethers: copy_examples, copy_bin, copy_files, install_examples, example_packages

pkg_dir = dirname(@__DIR__)

@testset "copy_files" begin
    mktempdir() do dir
        cd(dir) do
            files = copy_files("examples", ["Tether_01.jl"])
            @test files == ["Tether_01.jl"]
            @test isdir("examples")
            @test isfile(joinpath("examples", "Tether_01.jl"))
            @test read(joinpath(pkg_dir, "examples", "Tether_01.jl")) ==
                  read(joinpath("examples", "Tether_01.jl"))

            # overwrite=false must not touch an existing file
            dst = joinpath("examples", "Tether_01.jl")
            write(dst, "modified")
            copy_files("examples", ["Tether_01.jl"]; overwrite=false)
            @test read(dst, String) == "modified"

            # overwrite=true (the default) must replace it
            copy_files("examples", ["Tether_01.jl"])
            @test read(dst, String) != "modified"
        end
    end
end

@testset "copy_examples" begin
    mktempdir() do dir
        cd(dir) do
            copy_examples()
            @test isdir("examples")
            copied = readdir("examples")
            @test "Project.toml" ∉ copied
            @test !any(f -> startswith(f, "Manifest"), copied)
            src_files = filter(readdir(joinpath(pkg_dir, "examples"))) do file
                file != "Project.toml" && !startswith(file, "Manifest")
            end
            @test sort(copied) == sort(src_files)
        end
    end
end

@testset "example_packages" begin
    pkgs = example_packages()
    @test "Tethers" ∉ pkgs
    @test "LinearAlgebra" ∉ pkgs
    @test "ModelingToolkit" ∈ pkgs
    @test "GLMakie" ∈ pkgs
    @test length(pkgs) == length(unique(pkgs))
end

@testset "copy_bin" begin
    mktempdir() do dir
        cd(dir) do
            copy_bin()
            @test isdir("bin")
            # only `run_julia` is copied; the other scripts in `bin` are specific to a
            # clone of this repository
            @test readdir("bin") == ["run_julia"]
            @test read(joinpath(pkg_dir, "bin", "run_julia")) ==
                  read(joinpath("bin", "run_julia"))
        end
    end
end

@testset "install_examples" begin
    mktempdir() do dir
        cd(dir) do
            install_examples(false)
            @test isdir("examples")
            @test isdir("bin")
            copied = readdir("examples")
            @test "Project.toml" ∉ copied
            @test !any(f -> startswith(f, "Manifest"), copied)
            src_files = filter(readdir(joinpath(pkg_dir, "examples"))) do file
                file != "Project.toml" && !startswith(file, "Manifest")
            end
            @test sort(copied) == sort(src_files)
            @test readdir("bin") == ["run_julia"]
        end
    end
end
nothing
