function retrieve_surface(algorithm::DiscreteShape, img::AbstractArray, iterations::Int=2000)
    ρ,I,σ,τ = estimate_img_properties(img)
    E = Array{Float64}(img)
    E = E[1:2:end,1:2:end]
    M,N=size(E)
    p = zeros(axes(E))
    q = zeros(axes(E))
    δp = zeros(axes(E))
    δq = zeros(axes(E))
    R = zeros(axes(E))
    Z = zeros(axes(E))
    λ = 1000
    w = centered(0.25*[0 1 0;1 0 1;0 1 0])
    wx = zeros(axes(E))
    wy = zeros(axes(E))
    for i in CartesianIndices(wx)
        wx[i] = (2 * π * i[2]) / M
        wy[i] = (2 * π * i[1]) / N
    end
    return solve_EulerLagrange(ρ,I,iterations,δp,δq,w,p,q,R,λ,wx,wy,E,Z)
end
