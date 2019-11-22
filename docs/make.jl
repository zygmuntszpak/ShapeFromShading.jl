push!(LOAD_PATH,"../src/")
using Documenter
using ShapeFromShading

makedocs(
    sitename = "Documentation",
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Function Reference" => "reference.md"
    ]
    # modules = [ShapeFromShading]
)
deploydocs(repo = "github.com/betttris13/ShapeFromShading.jl.git")
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
