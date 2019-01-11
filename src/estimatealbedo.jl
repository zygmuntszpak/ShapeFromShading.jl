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
    τ = atan(mean(nEy)/mean(nEx))
    if τ < 0
        τ = τ + π;
    end
    I = [cos(τ)*sin(σ),sin(τ)*sin(σ),cos(σ)]
    return ρ,I,σ,τ
end
