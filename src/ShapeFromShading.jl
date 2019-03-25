module ShapeFromShading

using Images
using Statistics
using FFTW
using DSP
using LinearAlgebra

abstract type ShapeAlgorithm end
struct DiscreteShape <: ShapeAlgorithm end
struct DiscreteShapeBound <: ShapeAlgorithm end
struct Pentland <: ShapeAlgorithm end
struct Shah <: ShapeAlgorithm end
struct Photometric <: ShapeAlgorithm end

abstract type SynthShape end
struct SynthSphere <: SynthShape end
struct Ripple <: SynthShape end

include("common.jl")
include("syntheticsurface.jl")
include("estimatealbedo.jl")
include("discreteshape.jl")
include("discreteshapebound.jl")
include("pentland.jl")
include("shah.jl")
include("photometric.jl")

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
    Ripple
end # module
