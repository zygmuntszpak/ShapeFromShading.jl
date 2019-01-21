module ShapeFromShading

using Images
using Statistics
using AbstractFFTs
using DSP

abstract type ShapeAlgorithm end
struct DiscreteShape <: ShapeAlgorithm end
struct DiscreteShapeBound <: ShapeAlgorithm end
struct Pentland <: ShapeAlgorithm end
struct Shah <: ShapeAlgorithm end

include("common.jl")
include("syntheticsurface.jl")
include("estimatealbedo.jl")
include("discreteshape.jl")
include("discreteshapebound.jl")

export
    # main functions
    generate_surface,
    estimate_img_properties,
    retrieve_surface,
    DiscreteShape,
    DiscreteShapeBound
end # module
