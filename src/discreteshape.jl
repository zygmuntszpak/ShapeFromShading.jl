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
    for i = 1:iterations
        percent = (i/iterations)*100
        @show i,percent
        δp = imfilter(p,reflect(w),"replicate")
        δq = imfilter(q,reflect(w),"replicate")
        R = (ρ.*(-I[1] .* p .- I[2].* q .+ I[3]))./sqrt.(1 .+ p.^2 .+ q.^2)
        pq = (1 .+ p.^2 .+ q.^2)
        dR_dp = ((-ρ * I[1]) ./ sqrt.(pq)) .+ ((-I[1] * ρ) .* p .- I[2] * ρ .* q .+ I[3] * ρ) .* (-1 .* p .* (pq .^(-3/2)))
        dR_dq = ((-ρ*I[2]) ./ sqrt.(pq)) .+ (-I[1] * ρ .* p - I[2] * ρ .* q .+ I[3] * ρ) .* (-1 .* q .* (pq .^(-3/2)))
        p = δp .+ (1/(4*λ))*(E.-R).*dR_dp;
        q = δq .+ (1/(4*λ))*(E.-R).*dR_dq;
        Cp = fft(p)
        Cq = fft(q)
        C = -1im.*(wx .* Cp .+ wy .* Cq)./(wx.^2 .+ wy.^2)
        Z = abs.(ifft(C))
        p = ifft(1im .* wx .* C);
        q = ifft(1im .* wy .* C);
    end
    return Z
end
