using Tethers, Pkg
if ! ("Documenter" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using Documenter: deploydocs, makedocs

src="src"
dst="docs/src/src"
mkpath(dst)
cp(src, dst; force=true)
src="docs/images"
dst="docs/src/docs/images"
mkpath(dst)
cp(src, dst; force=true)

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
        "Theory" => "theory.md",
        "Examples" => "examples.md",
        "VSCode IDE" => "vscode.md",
        "Python and Julia" => "python.md",
        "References" => "references.md"
    ])
deploydocs(repo = "github.com/ufechner7/Tethers.jl.git")
