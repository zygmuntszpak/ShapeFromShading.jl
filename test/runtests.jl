using Images, ShapeFromShading
using Test

@testset "ShapeFromShading.jl" begin
    include("syntheticsurface.jl")
    include("shah.jl")
    include("pentland.jl")
    include("discreteshape.jl")
    include("discreteshapebound.jl")
    include("estimatealbedo.jl")
end
