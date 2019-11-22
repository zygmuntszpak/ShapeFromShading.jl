push!(LOAD_PATH,"../src/")
using Documenter
using ShapeFromShading

makedocs(
    sitename = "Documentation",
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md"
    ],
    modules = [ShapeFromShading]
)
deploydocs(repo = "github.com/betttris13/ShapeFromShading.jl.git")
