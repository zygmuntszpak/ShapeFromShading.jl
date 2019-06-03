using Images, ShapeFromShading
using Test

@testset "ShapeFromShading" begin
    include("shah.jl")
    include("pentland.jl")
    include("discreteshape.jl")
    include("discreteshapebound.jl")
    include("photometric.jl")
end

@testset "Integration" begin
    include("frankot.jl")
    include("path.jl")
    include("splitpath.jl")
    include("horn.jl")
    include("durou.jl")
end

@testset "SyntheticSurface" begin
    include("synthsphere.jl")
    include("synthvase.jl")
    include("tent.jl")
    include("dem.jl")
    include("cake.jl")
    include("cake2.jl")
    include("ripple.jl")
    include("ripple2.jl")
    include("cup.jl")
    include("synthgaussian.jl")
    include("supergaussian.jl")
    include("estimatealbedo.jl")
end
