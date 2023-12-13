using Tethers, Pkg
if ("TestEnv" ∈ keys(Pkg.project().dependencies))
    if ! ("Documenter" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
end
using Documenter: deploydocs, makedocs

src="README.md"
dst="docs/src/index.md"
if ! isfile(dst)
    cp(src, dst)
elseif readlines(src) != readlines(dst)
    cp(src, dst; force=true)
end

makedocs(;
    authors="Uwe Fechner <fechner@aenarete.eu>",
    sitename = "Tethers.jl", 
    modules = [Tethers], 
    doctest = false,
    pages=[
        "Readme" => "index.md",
        "Python and Julia" => "python.md"
    ])
deploydocs(repo = "github.com/ufechner7/PkgHelpers.jl.git")
