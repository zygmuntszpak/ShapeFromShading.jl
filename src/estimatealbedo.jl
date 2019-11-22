@doc raw"""
```
albedo,illumination_direction,slant,tilt = estimate_img_properties(img::AbstractArray)
```
Attempts to calculate the illumination properties of an image for use in SFS
algorithms.
# Output
Returns albedo (ρ), illumination vector (I), slant (σ) and tilt (τ) of the image as Float,
Array{Float} of size 3, Float and Float respectively.
# Details
The illumination properties can be estimated using the first and second
moments and the gradients of the image (``E_x, E_y``) brightness as per the following:

```math
\begin{gathered}
\mu_1=\dfrac{\pi}{4}\rho\cos\sigma\\\\\mu_2=\dfrac{1}{6}\rho^2(1+3\cos^2\sigma)\\
\\\rho=\dfrac{\sqrt{6\pi\mu_2-48\mu_1^2}}{\pi}\\\\\cos\sigma=\dfrac{4\mu_1}{\sqrt
{6\pi\mu_2-48\mu_1^2}}\\\\\tan\tau=\dfrac{\text{average}(E_y)}{\text{average}
(E_x)}\\\\I=[\cos\tau\sin\sigma,\sin\tau\sin\sigma,\cos\sigma]
\end{gathered}
```
Note: Properties are an estimate and may not be correct under all conditions.
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
img = generate_surface(SynthSphere(), 1, [0.2,0,0.9], radius = 5)

#estimate the properties
albedo, illumination_direction, slant, tilt = estimate_img_properties(img)
```
# Reference
1. T. Ping-Sing and M. Shah, "Shape from shading using linear approximation", Image and Vision Computing, vol. 12, no. 8, pp. 487-498, 1994. [doi:10.1016/0262-8856(94)90002-7](https://doi.org/10.1016/0262-8856(94)90002-7)
"""
function estimate_img_properties(img::AbstractArray)
    E = Array{Float64}(img)
    #calculate spatial gradient of img
    Ey, Ex = imgradients(E, KernelFactors.sobel, "replicate")
    #normalize gradient
    nEx = similar(E)
    nEy = similar(E)
    Exy = sqrt.((Ex.^2) .+ (Ey.^2))
    nEx = Ex ./ (Exy .+ eps())
    nEy = Ey ./ (Exy .+ eps())

    #calculate means
    μ₁ = mean(E)
    μ₂ = mean(E.^2)

    #calculate desired values
    g = sqrt((6.0 * (π^2) * μ₂) - (48.0 * (μ₁^2)))
    ρ = g / π
    t = (4 * μ₁) / g
    if t > 1.0
        t = 1.0
    elseif t < -1.0
        t = -1.0
    end
    σ = acos(t)
    τ = atan(mean(nEy) / mean(nEx))
    if τ < 0
        τ = τ + π;
    end
    I = [cos(τ) * sin(σ), sin(τ) * sin(σ), cos(σ)]
    return ρ, I, σ, τ
end
