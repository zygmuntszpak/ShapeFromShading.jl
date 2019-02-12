@doc raw"""
```
Z,p,q = retrieve_surface(algorithm::Photometric, img1::AbstractArray, img2::AbstractArray, img3::AbstractArray, illumination_direction1::Vector{T} where T <: Real, illumination_direction2::Vector{T} where T <: Real, illumination_direction3::Vector{T} where T <: Real; Integration_Scheme::IntegrationScheme = Horn())
```
Attempts to produce the gradients and integrate them to retiereve a height map
using three-point photometric stereo.

This algorithm solves a system of linear equation formed by varying the lighting
conditions on an object while holding the viewing angle constant. Note that
the ilumination direction id defined negative of most other algorithms and as
such inouts may need to be corrected (see example).
# Output
Returns an M by N array (matching dimensions of original image) of Float `Z`
that represents the reconstructed height at the point and the gradients in
M by N arrays of Float `p` and `q`.
# Details
Let ⁠ñ₁, ñ₂ and ñ₃ define the ilumination direction let N be a 3 X 3 matrix
whos rows are formed by the row vectors of the illumination direction as per:
```math
N=\begin{bmatrix}&n_{11}&n_{12}&n_{13}\\&n_{21}&n_{22}&n_{23}\\&n_{31}&n_{32}&n_{33}\end{bmatrix}
```
and I₁, I₂ and I₃ be the intensity values at a point (x,y) and let Ĩ=[I₁,I₂,I₃]′
be the column vector formed by these values.

The surface normal at the point (x,y) can now be directly calculated using:
```math
ñ=N^{-1}Ĩ
```
where n is the surface normal.
# Arguments
The function arguments are described in more detail below.
##  `img1`
An `AbstractArray` storing the grayscale value of each pixel within
the range [0,1].
##  `img2
An `AbstractArray` storing the grayscale value of each pixel within
the range [0,1].
##  `img3`
An `AbstractArray` storing the grayscale value of each pixel within
the range [0,1].
## `illumination_direction1`
A `Vector{T} where T <: Real` that specifies the first illumination vector to be used by the
algorithm. The `illumination_direction` should be a vector of the form [x,y,z]
where x,y,z are int he range [0,1].
## `illumination_direction2`
A `Vector{T} where T <: Real` that specifies the second illumination vector to be used by the
algorithm. The `illumination_direction` should be a vector of the form [x,y,z]
where x,y,z are int he range [0,1].
## `illumination_direction3`
A `Vector{T} where T <: Real` that specifies the third illumination vector to be used by the
algorithm. The `illumination_direction` should be a vector of the form [x,y,z]
where x,y,z are int he range [0,1].
## `Integration_Scheme`
A keywork argument of type `IntegrationScheme` which specifies which integration
scheme is to be used to reconstruct the surface from the calculated normals.
```julia
using Images, Makie, ShapeFromShading

#generate synthetic images
n₁ = [0,0,1.0]
n₂ = [1.0,0.2,1.0]
n₃ = [0.2,0.9,1.0]
img1, img2, img3 = generate_photometric(n₁, n₂, n₃, 1, shape = SynthSphere(), r = 5)

#calculate the heightmap
Z,p,q = retrieve_surface(Photometric(), img1, img2, img3, -n₁, -n₂, -n₃, Integration_Scheme = Durou())

#normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(Z)

#display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:2
surface(r, r, Z)
```
# Reference
1. R. Woodham, "Photometric Method For Determining Surface Orientation From Multiple Images", Optical Engineering, vol. 19, no. 1, 1980. [doi:10.1117/12.7972479](https://doi.org/10.1117/12.7972479)
"""
function retrieve_surface(algorithm::Photometric, img1::AbstractArray, img2::AbstractArray, img3::AbstractArray, illumination_direction1::Vector{T} where T <: Real, illumination_direction2::Vector{T} where T <: Real, illumination_direction3::Vector{T} where T <: Real; Integration_Scheme::IntegrationScheme = Horn())
    #setup illumination matrix
    N = zeros(Float64, 3, 3)
    for i = 1:3
        N[1,i] = illumination_direction1[i]
        N[2,i] = illumination_direction2[i]
        N[3,i] = illumination_direction3[i]
    end
    Ninv = inv(N)

    #setup image and gradient arrays
    I₁ = Float64.(img1)
    I₂ = Float64.(img2)
    I₃ = Float64.(img3)
    p = zeros(Float64, size(I₁))
    q = zeros(Float64, size(I₁))

    #calculate gradients
    for i in CartesianIndices(I₁)
        I = [I₁[i], I₂[i], I₃[i]]
        n = Ninv * I
        if n[3] == 0.0
            n = [0.0,0.0,-1.0]
        end
        p[i] = -n[1] / n[3]
        q[i] = -n[2] / n[3]
    end

    #reconstruct surface
    Z = convert_gradient(Integration_Scheme, p, q)
    return Z, p, q
end
