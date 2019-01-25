"""
```
albedo,illumination_direction,slant,tilt = estimate_img_properties(img::AbstractArray)
```
Attempts to calculate the propeties of an image for use in shape form shading
algarithms.
# Output
Returns albedo , ilumination vector, slant and tilt of the image as Float,
Float, Array{Float} of size 3 and Float respectivly.
# Details
propeties are an estimate and may not be correct under all conditions.
# Arguments
The function arguments are described in more detail below.
##  `img`
An `AbstractArray` storing the grayscale value of each pixel within
the range [0,1].
# Example
Compute properties of synthetic image generated using `generate_surface`.
```julia
using Images, ShapeFromShading

#generate synthetic image
img = generate_surface(0.5, [0.2,0,0.9], radius = 5)

#estimate the properties
albedo,illumination_direction,slant,tilt = estimate_img_properties(img)
```
# Reference
1. T. Ping-Sing and M. Shah, "Shape from shading using linear approximation", Image and Vision Computing, vol. 12, no. 8, pp. 487-498, 1994. [doi:10.1016/0262-8856(94)90002-7](https://doi.org/10.1016/0262-8856(94)90002-7)
"""
function estimate_img_properties(img::AbstractArray)
    E = Array{Float64}(img)
    #calculate spatial gradient of img
    Ey,Ex = imgradients(E, KernelFactors.sobel, "replicate")
    #normalize gradient
    nEx = similar(E)
    nEy = similar(E)
    Exy = sqrt.((Ex.^2) .+ (Ey.^2))
    nEx = Ex./(Exy.+eps())
    nEy = Ey./(Exy.+eps())

    #calculate means
    μ₁ = mean(E)
    μ₂ = mean(E.^2)

    #calculate desired values
    g = sqrt(6*(π^2)*μ₂-48*(μ₁^2))
    ρ = g/π
    σ = acos((4*μ₁)/g)
    #@show mean(nEy),mean(nEx)
    τ = atan(mean(nEy)/mean(nEx))
    if τ < 0
        τ = τ + π;
    end
    I = [cos(τ)*sin(σ),sin(τ)*sin(σ),cos(σ)]
    @show typeof(I)
    return ρ,I,σ,τ
end
