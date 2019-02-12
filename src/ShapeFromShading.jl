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

abstract type ShapeAlgorithm end
struct Determanistic <: ShapeAlgorithm end
struct Determanistic2 <: ShapeAlgorithm end
struct DiscreteShape <: ShapeAlgorithm end
struct DiscreteShapeBound <: ShapeAlgorithm end
struct Pentland <: ShapeAlgorithm end
struct Photometric <: ShapeAlgorithm end
struct Shah <: ShapeAlgorithm end

abstract type IntegrationScheme end
struct Durou <: IntegrationScheme end
struct Frankot <: IntegrationScheme end
struct Horn <: IntegrationScheme end
struct Path <: IntegrationScheme end
struct SplitPath <: IntegrationScheme end


abstract type SynthShape end
struct Cake <: SynthShape end
struct Cake2 <: SynthShape end
struct Cup <: SynthShape end
struct Dem <: SynthShape end
struct Ripple <: SynthShape end
struct Ripple2 <: SynthShape end
struct SynthGaussian <: SynthShape end
struct SuperGaussian <: SynthShape end
struct SynthSphere <: SynthShape end
struct SynthVase <: SynthShape end
struct Tent <: SynthShape end

include("benchmark.jl")
include("common.jl")
include("estimatealbedo.jl")
include("discreteshape.jl")
include("discreteshapebound.jl")
include("determanistic.jl")
include("integration.jl")
include("pentland.jl")
include("photometric.jl")
include("shah.jl")
include("syntheticsurface.jl")

export
    # main functions and datatypes
    benchmark_iterative,
    benchmark_noniterative,
    benchmark_shape,
    Cake,
    Cake2,
    compare_benchmark,
    convert_gradient,
    createImage,
    Cup,
    Dem,
    Determanistic,
    Determanistic2,
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
    SplitPath,
    SuperGaussian,
    sythetic_gradient,
    SynthGaussian,
    SynthSphere,
    SynthVase,
    Tent
end # module
