function retrieve_surface(algorithm::DiscreteShapeBound, img::AbstractArray, iterations::Int=2000)
    ρ,I,σ,τ = estimate_img_properties(img)
    return retrieve_surface(DiscreteShapeBound(), img, ρ, I, iterations)
end

function retrieve_surface(algorithm::DiscreteShapeBound, img::AbstractArray, albedo::Real, illumination_direction::Vector{T} where T <: Real, iterations::Int=2000)
    ρ = albedo
    I = illumination_direction
    E = Array{Float64}(img)
    E = E[1:2:end,1:2:end]
    M,N=size(E)
    p = zeros(axes(E))
    q = zeros(axes(E))
    δp = zeros(axes(E))
    δq = zeros(axes(E))
    R = zeros(axes(E))
    Z = zeros(axes(E))
    for i in CartesianIndices(E)
        if E[i]>0.75
            Z[i]=-100*E[i]
        else
            Z[i]=0
        end
    end
    q,p = imgradients(Z, KernelFactors.sobel, "replicate")
    λ = 1000
    w = centered(0.25*[0 1 0;1 0 1;0 1 0])
    wx,wy = setup_transform_values(M,N)
    return solve_EulerLagrange(ρ,I,iterations,δp,δq,w,p,q,R,λ,wx,wy,E,Z)
end
