# ShapeFromShading.jl

The following documents the functionality currently available in ShapeFormShading.jl.

## Synthetic data generation:

```@docs
generate_surface(::AbstractSyntheticShape, ::Real, ::Vector{T} where T <: Real; ::Int, ::Real, ::Bool)
```

## Normal Integration:

```@docs
AbstractIntegrator
```

```@docs
Frankot(::AbstractArray, ::AbstractArray)
```

```@docs
Path(::AbstractArray, ::AbstractArray)
```

```@docs
SplitPath(::AbstractArray, ::AbstractArray)
```

```@docs
Horn(::AbstractArray, ::AbstractArray)
```

```@docs
Durou(::AbstractArray, ::AbstractArray)
```

```@docs
Quadratic(::AbstractArray, ::AbstractArray)
```

```@docs
TotalVariation(::AbstractArray, ::AbstractArray)
```

```@docs
NonConvex1(::AbstractArray, ::AbstractArray)
```

```@docs
NonConvex2(::AbstractArray, ::AbstractArray)
```

```@docs
AnisotropicDiffusion(::AbstractArray, ::AbstractArray)
```

```@docs
MumfordShah(::AbstractArray, ::AbstractArray)
```

## Shape From Shading:

```@docs
DiscreteShape(::AbstractArray)
```

```@docs
DiscreteShapeBound(::AbstractArray)
```

```@docs
Pentland(::AbstractArray)
```

```@docs
Photometric(::AbstractArray,::AbstractArray,::AbstractArray)
```

```@docs
Shah(::AbstractArray)
```

```@docs
Deterministic(::AbstractArray)
```

```@docs
DeterministicCustom(::AbstractArray)
```

## Benchmarking:



## Miscellaneous:
