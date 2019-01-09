module ShapeFromShading

using Images

abstract type ShapeAlgorithm end
struct DiscreteShape <: ShapeAlgorithm end
struct DiscreteShapeBound <: ShapeAlgorithm end
struct Pentland <: ShapeAlgorithm end
struct Shah <: ShapeAlgorithm end

include("discreteshape.jl")
include("syntheticsurface.jl")

export
    # main functions
    generate_surface,
    retrieve_surface,
    DiscreteShape
end # module
