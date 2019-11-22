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
@with_kw struct Cake <: AbstractSyntheticShape
    radius::Real = 50
end
@with_kw struct Cake2 <: AbstractSyntheticShape
    radius::Real = 75
end
@with_kw struct Dem <: AbstractSyntheticShape
    radius::Real = 75
end
@with_kw struct Prism <: AbstractSyntheticShape
    radius::Real = 75
end
@with_kw struct Ripple <: AbstractSyntheticShape
    radius::Real = 10
end
@with_kw struct Ripple2 <: AbstractSyntheticShape
    radius::Real = 7.5
end
@with_kw struct SynthGaussian <: AbstractSyntheticShape
    radius::Real = 75
end
@with_kw struct SuperGaussian <: AbstractSyntheticShape
    radius::Real = 100
end
@with_kw struct SynthSphere <: AbstractSyntheticShape
    radius::Real = 50
end
@with_kw struct SynthVase <: AbstractSyntheticShape
    radius::Real = 78.5
end
@with_kw struct Tent <: AbstractSyntheticShape
    radius::Real = 50
end

@doc raw"""
```
(scheme::AbstractIntegrator)(p::AbstractArray, q::AbstractArray)
```
Integrates the gradient field defined by `p` and `q` using various methods as
described in the following sections.
# Output
Returns a heightmap `Z` which defines the reconstructed surface.
# Details
Using the method defined in each type of integrator, each algorithm is run using
parameters passed to the struct when creating the integrator. Results can vary
wildly depending on the chosen integrator.
Available methods are; Frankot, Path, SplitPath, Horn, Durou, Quadratic,
TotalVariation, NonConvex1, NonConvex2, AnisotropicDiffusion and MumfordShah.
See individual documentation for each method for more detials.
# Parameters
Each algorithm has its own set of parameters which are discussed in each the section
for each method.
# Arguments
## `p` and `q`:
`AbstractArray`'s which define the gradient field to be integrated.
# Example
The following example demonstrates  the use of the `Frankot` integrator.
```julia
using ShapeFromShading, Makie

# Generate synthetic gradients
p, q = synthetic_gradient(SynthSphere(38), img_size = 151)

# Create a Frankot() integrator
frankot = Frankot()

# Calculate the heightmap from the gradients
Z = frankot(p, q)

# Normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(Z)

# Display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
surface(r, r, Z)
```
"""
abstract type AbstractIntegrator end
@with_kw struct AnisotropicDiffusion <: AbstractIntegrator
    z::AbstractArray
    μ::Real = 5.0
    ν::Real = 10.0
    λ::AbstractArray = fill(10.0^-6, size(z))
    mask::AbstractArray = fill(1.0, size(z))
    max_iter::Int = 10
end
@with_kw struct Durou <: AbstractIntegrator
    ϵ::Real = 1.0
    max_iter::Real = 1000
end
struct Frankot <: AbstractIntegrator end
@with_kw struct Horn <: AbstractIntegrator
    ϵ::Real = 1.0
    max_iter::Real = 10000
end
@with_kw struct MumfordShah <: AbstractIntegrator
    z::AbstractArray
    μ::Real = 10.0
    ϵ::Real = 0.1
    λ::AbstractArray = fill(10.0^-6, size(z))
    mask::AbstractArray = fill(1.0, size(z))
    max_iter::Int = 50
end
@with_kw struct NonConvex1 <: AbstractIntegrator
    z::AbstractArray
    β::Real = 0.5
    λ::AbstractArray = fill(10.0^-6, size(z))
    mask::AbstractArray = fill(1.0, size(z))
    max_iter::Int = 100
end
@with_kw struct NonConvex2 <: AbstractIntegrator
    z::AbstractArray
    γ::Real = 1.0
    λ::AbstractArray = fill(10.0^-6, size(z))
    mask::AbstractArray = fill(1.0, size(z))
    max_iter::Int = 100
end
struct Path <: AbstractIntegrator end
struct SplitPath <: AbstractIntegrator end
@with_kw struct TotalVariation <: AbstractIntegrator
    z::AbstractArray
    α::Real = 1.0
    λ::AbstractArray = fill(10.0^-6, size(z))
    mask::AbstractArray = fill(1.0, size(z))
    max_iter::Int = 100
end
@with_kw struct Quadratic <: AbstractIntegrator
    z::AbstractArray
    λ::AbstractArray = fill(10.0^-6, size(z))
    mask::AbstractArray = fill(1.0, size(z))
end

abstract type AbstractShapeAlgorithm end
@with_kw struct Deterministic <: AbstractShapeAlgorithm
     p::AbstractArray
     q::AbstractArray
     integration_factor::Real = 10
     smoothness::Real = 20
     β::Real = 1.0
     integration_scheme::AbstractIntegrator = Horn()
end
@with_kw struct DeterministicCustom <: AbstractShapeAlgorithm
    p::AbstractArray
    q::AbstractArray
    integration_factor::Real = 10
    smoothness::Real = 20
    β::Real = 1.0
    integration_scheme::AbstractIntegrator = Horn()
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
    integration_scheme::AbstractIntegrator = Horn()
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
    AbstractIntegrator,
    AbstractShapeAlgorithm,
    AbstractSyntheticShape,
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
    NonConvex2,
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
