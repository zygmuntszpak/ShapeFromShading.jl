module ShapeFromShading

using AlgebraicMultigrid
using Blink
using Distributions
using DSP
using FFTW
using Images
using Interact
using IterativeSolvers
using LinearAlgebra
using Makie
using Optim
using Parameters
using Preconditioners
using SparseArrays
using Statistics

abstract type AbstractSyntheticShape end
struct Cake <: AbstractSyntheticShape end
struct Cake2 <: AbstractSyntheticShape end
struct Cup <: AbstractSyntheticShape end
struct Dem <: AbstractSyntheticShape end
struct Prism <: AbstractSyntheticShape end
struct Ripple <: AbstractSyntheticShape end
struct Ripple2 <: AbstractSyntheticShape end
struct SynthGaussian <: AbstractSyntheticShape end
struct SuperGaussian <: AbstractSyntheticShape end
struct SynthSphere <: AbstractSyntheticShape end
struct SynthVase <: AbstractSyntheticShape end
struct Tent <: AbstractSyntheticShape end

abstract type AbstractIntegrationScheme end  # Rename to AbstractIntegrator?
@with_kw struct AnisotropicDiffusion <: AbstractIntegrationScheme
    z::AbstractArray
    μ::Real = 5.0
    ν::Real = 10.0
    λ::AbstractArray = fill(10.0^-6, size(z))
    mask::AbstractArray = fill(1.0, size(z))
    max_iter::Int = 10
end
@with_kw struct Durou <: AbstractIntegrationScheme
    ϵ::Real = 1.0
    max_iter::Real = 1000
end
struct Frankot <: AbstractIntegrationScheme end
@with_kw struct Horn <: AbstractIntegrationScheme
    ϵ::Real = 1.0
    max_iter::Real = 10000
end
@with_kw struct MumfordShah <: AbstractIntegrationScheme
    z::AbstractArray
    μ::Real = 10.0
    ϵ::Real = 0.1
    λ::AbstractArray = fill(10.0^-6, size(z))
    mask::AbstractArray = fill(1.0, size(z))
    max_iter::Int = 50
end
@with_kw struct NonConvex1 <: AbstractIntegrationScheme
    z::AbstractArray
    β::Real = 0.5
    λ::AbstractArray = fill(10.0^-6, size(z))
    mask::AbstractArray = fill(1.0, size(z))
    max_iter::Int = 100
end
struct Path <: AbstractIntegrationScheme end
struct SplitPath <: AbstractIntegrationScheme end
@with_kw struct TotalVariation <: AbstractIntegrationScheme
    z::AbstractArray
    α::Real = 1.0
    λ::AbstractArray = fill(10.0^-6, size(z))
    mask::AbstractArray = fill(1.0, size(z))
    max_iter::Int = 100
end
@with_kw struct Quadratic <: AbstractIntegrationScheme
    z::AbstractArray
    λ::AbstractArray = fill(10.0^-6, size(z))
    mask::AbstractArray = fill(1.0, size(z))
end

abstract type AbstractShapeAlgorithm end
@with_kw struct Deterministic <: AbstractShapeAlgorithm
     pin::AbstractArray
     qin::AbstractArray
     integration_factor::Real = 10
     smoothness::Real = 20
end
@with_kw struct DeterministicCustom <: AbstractShapeAlgorithm
     pin::AbstractArray
     qin::AbstractArray
     integration_factor::Real = 10
     smoothness::Real = 20
end
@with_kw struct DiscreteShape <: AbstractShapeAlgorithm
    iterations::Int = 2000
    smoothness::Int = 1000
    albedo::Real = Inf
    illumination_direction::Vector{T} where T <: Real = [Inf, Inf, Inf]
end
@with_kw struct DiscreteShapeBound <: AbstractShapeAlgorithm
    iterations::Int = 2000
    smoothness::Int = 1000
    albedo::Real = Inf
    illumination_direction::Vector{T} where T <: Real = [Inf, Inf, Inf]
end
@with_kw struct MultiResolutionHybrid <: AbstractShapeAlgorithm
    T₀::Real = 500.0
    α::Real = 0.99
    integration_factor::Real = 10
    smoothness::Real = 20
    integration_factor2::Real = 10
    smoothness2::Real = 20
end
@with_kw struct Pentland <: AbstractShapeAlgorithm
    slant::Real = Inf
    tilt::Real = Inf
end
@with_kw struct Photometric <: AbstractShapeAlgorithm
    illumination_direction1::Vector{T} where T <: Real = [Inf, Inf, Inf]
    illumination_direction2::Vector{T} where T <: Real = [Inf, Inf, Inf]
    illumination_direction3::Vector{T} where T <: Real = [Inf, Inf, Inf]
    integration_scheme::AbstractIntegrationScheme = Horn()
end
@with_kw struct Shah <: AbstractShapeAlgorithm
    slant::Real = Inf
    tilt::Real = Inf
    iterations::Int = 200
end
@with_kw struct SimulatedAnnealing <: AbstractShapeAlgorithm
    pin::AbstractArray
    qin::AbstractArray
    T₀::Real = 500.0
    α::Real = 0.99
    integration_factor::Real = 10
    smoothness::Real = 20
    max_iter::Int = 40000
    debug::Bool = false
end

include("benchmark.jl")
include("common.jl")
include("estimatealbedo.jl")
include("discreteshape.jl")
include("discreteshapebound.jl")
include("deterministic.jl")
include("integration.jl")
include("multiresolutionhybrid.jl")
include("pentland.jl")
include("photometric.jl")
include("shah.jl")
include("syntheticsurface.jl")
include("simulatedannealing.jl")

export
    # main functions and datatypes
    AnisotropicDiffusion,
    benchmark_image,
    benchmark_integration,
    benchmark_normals,
    benchmark_shape,
    Cake,
    Cake2,
    convert_gradient,
    create_image,
    Cup,
    Dem,
    Deterministic,
    DeterministicCustom,
    DiscreteShape,
    DiscreteShapeBound,
    Durou,
    estimate_img_properties,
    Frankot,
    gen_mask,
    generate_photometric,
    generate_surface,
    generate_surface_data,
    ground_truth,
    Horn,
    MultiResolutionHybrid,
    MumfordShah,
    NonConvex1,
    Path,
    Pentland,
    Photometric,
    Prism,
    retrieve_surface,
    Ripple,
    Ripple2,
    Shah,
    SimulatedAnnealing,
    SplitPath,
    SuperGaussian,
    synthetic_gradient,
    SynthGaussian,
    SynthSphere,
    SynthVase,
    Tent,
    TotalVariation,
    Quadratic,
    MumfordShah2
end # module
