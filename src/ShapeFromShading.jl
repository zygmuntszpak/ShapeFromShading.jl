module ShapeFromShading

using Blink
using Distributions
using DSP
using FFTW
using Images
using Interact
using LinearAlgebra
using Makie
using Optim
using Statistics

abstract type AbstractShapeAlgorithm end
struct Deterministic <: AbstractShapeAlgorithm end
struct DeterministicCustom <: AbstractShapeAlgorithm end
struct DiscreteShape <: AbstractShapeAlgorithm end
struct DiscreteShapeBound <: AbstractShapeAlgorithm end
struct Pentland <: AbstractShapeAlgorithm end
struct Photometric <: AbstractShapeAlgorithm end
struct Shah <: AbstractShapeAlgorithm end
struct SimulatedAnnealing <: AbstractShapeAlgorithm end

abstract type AbstractIntegrationScheme end
struct Durou <: AbstractIntegrationScheme end
struct Frankot <: AbstractIntegrationScheme end
struct Horn <: AbstractIntegrationScheme end
struct Path <: AbstractIntegrationScheme end
struct SplitPath <: AbstractIntegrationScheme end


abstract type AbstractSyntheticShape end
struct Cake <: AbstractSyntheticShape end
struct Cake2 <: AbstractSyntheticShape end
struct Cup <: AbstractSyntheticShape end
struct Dem <: AbstractSyntheticShape end
struct Ripple <: AbstractSyntheticShape end
struct Ripple2 <: AbstractSyntheticShape end
struct SynthGaussian <: AbstractSyntheticShape end
struct SuperGaussian <: AbstractSyntheticShape end
struct SynthSphere <: AbstractSyntheticShape end
struct SynthVase <: AbstractSyntheticShape end
struct Tent <: AbstractSyntheticShape end

include("benchmark.jl")
include("common.jl")
include("estimatealbedo.jl")
include("discreteshape.jl")
include("discreteshapebound.jl")
include("deterministic.jl")
include("integration.jl")
include("pentland.jl")
include("photometric.jl")
include("shah.jl")
include("syntheticsurface.jl")
include("simulatedannealing.jl")

export
    # main functions and datatypes
    benchmark_iterative,
    benchmark_noniterative,
    benchmark_shape,
    Cake,
    Cake2,
    compare_benchmark,
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
    generate_photometric,
    generate_surface,
    generate_surface_data,
    ground_truth,
    Horn,
    Path,
    Pentland,
    Photometric,
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
    Tent
end # module
