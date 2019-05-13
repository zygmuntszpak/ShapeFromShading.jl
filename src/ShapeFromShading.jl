module ShapeFromShading

using Images
using Statistics
using FFTW
using DSP
using LinearAlgebra
using Distributions
using Makie
using Interact
using Blink

abstract type ShapeAlgorithm end
struct DiscreteShape <: ShapeAlgorithm end
struct DiscreteShapeBound <: ShapeAlgorithm end
struct Pentland <: ShapeAlgorithm end
struct Shah <: ShapeAlgorithm end
struct Photometric <: ShapeAlgorithm end

abstract type IntegrationScheme end
struct Frankot <: IntegrationScheme end
struct Path <: IntegrationScheme end
struct SplitPath <: IntegrationScheme end
struct Horn <: IntegrationScheme end
struct Durou <: IntegrationScheme end

abstract type SynthShape end
struct SynthSphere <: SynthShape end
struct Ripple <: SynthShape end
struct Ripple2 <: SynthShape end
struct Cake <: SynthShape end
struct Cake2 <: SynthShape end
struct Cup <: SynthShape end
struct SynthVase <: SynthShape end
struct Tent <: SynthShape end
struct Dem <: SynthShape end
struct SynthGaussian <: SynthShape end
struct SuperGaussian <: SynthShape end

include("common.jl")
include("syntheticsurface.jl")
include("estimatealbedo.jl")
include("discreteshape.jl")
include("discreteshapebound.jl")
include("pentland.jl")
include("shah.jl")
include("photometric.jl")
include("integration.jl")
include("benchmark.jl")

export
    # main functions
    generate_surface,
    estimate_img_properties,
    retrieve_surface,
    DiscreteShape,
    DiscreteShapeBound,
    Pentland,
    Shah,
    Photometric,
    convert_gradient,
    generate_photometric,
    sythetic_gradient,
    SynthSphere,
    Ripple,
    Ripple2,
    Cake,
    Cake2,
    ground_truth,
    Cup,
    SynthVase,
    Tent,
    Dem,
    Frankot,
    Path,
    SplitPath,
    Horn,
    Durou,
    benchmark_iterative,
    benchmark_noniterative,
    compare_benchmark,
    SynthGaussian,
    SuperGaussian,
    createImage
end # module
