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

@testset "Integration A-N" begin
    include("anisotropicdiffusion.jl")
    include("durou.jl")
    include("frankot.jl")
    include("horn.jl")
    include("mumfordshah.jl")
    nclude("nonconvex1.jl")
end

@testset "Integration O-Z" begin
    include("path.jl")
    include("splitpath.jl")
    include("totalvariation.jl")
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
