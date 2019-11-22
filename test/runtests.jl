using Images, ShapeFromShading
using Test

@testset "ShapeFromShading" begin
    include("deterministic.jl")
    include("deterministic_custom.jl")
    include("discreteshape.jl")
    include("discreteshapebound.jl")
    include("multiresolutionhybrid.jl")
    include("pentland.jl")
    include("photometric.jl")
    include("shah.jl")
    # include("SimulatedAnnealing.jl")
end

@testset "Anisotropicdiffusion" begin
    include("anisotropicdiffusion.jl")
end
@testset "Durou" begin
    include("durou.jl")
end
@testset "Frankot" begin
    include("frankot.jl")
end
@testset "Horn" begin
    include("horn.jl")
end
@testset "Mumfordshah" begin
    include("mumfordshah.jl")
end
@testset "Nonconvex1" begin
    include("nonconvex1.jl")
end
@testset "Path" begin
    include("path.jl")
end
@testset "Splitpath" begin
    include("splitpath.jl")
end
@testset "Totalvariation" begin
    include("totalvariation.jl")
end
@testset "Quadratic" begin
    include("quadratic.jl")
end

@testset "SyntheticSurface" begin
    include("cake.jl")
    include("cake2.jl")
    include("dem.jl")
    include("prism.jl")
    include("ripple.jl")
    include("ripple2.jl")
    include("supergaussian.jl")
    include("synthsphere.jl")
    include("synthvase.jl")
    include("synthgaussian.jl")
    include("tent.jl")
end

@testset "ImageGeneration" begin
    include("estimatealbedo.jl")
    include("generate_surface.jl")
end

@testset "OtherFunctions" begin
    include("mask.jl")
end

@testset "Benchmarks" begin
    include("benchmark_shape.jl")
    include("benchmark_integration.jl")
    include("benchmark_normals.jl")
    include("benchmark_image.jl")
end
