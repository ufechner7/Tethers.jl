using Test
using Tethers: copy_examples, copy_bin, copy_files, install_examples

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
            src_files = readdir(joinpath(pkg_dir, "examples"))
            @test sort(readdir("examples")) == sort(src_files)
        end
    end
end

@testset "copy_bin" begin
    mktempdir() do dir
        cd(dir) do
            copy_bin()
            @test isdir("bin")
            copied = readdir("bin")
            @test !any(f -> endswith(f, ".so"), copied)
            src_files = filter(f -> !endswith(f, ".so"), readdir(joinpath(pkg_dir, "bin")))
            @test sort(copied) == sort(src_files)
        end
    end
end

@testset "install_examples" begin
    mktempdir() do dir
        cd(dir) do
            install_examples(false)
            @test isdir("examples")
            @test isdir("bin")
            @test sort(readdir("examples")) == sort(readdir(joinpath(pkg_dir, "examples")))
            @test !any(f -> endswith(f, ".so"), readdir("bin"))
        end
    end
end
nothing
