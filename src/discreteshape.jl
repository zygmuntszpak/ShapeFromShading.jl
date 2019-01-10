function retrieve_surface(algorithm::DiscreteShape, img::AbstractArray)
    E = Array{Float64}(img)
    E = E[1:2:end,1:2:end]
    #ρ,I,σ,τ = estimate_img_properties(img)
    ρ,I = 0.884985933916937, [0.597336100966117,0.728954641551525,0.334387070687676]
    M,N=size(E)
    p = similar(E)
    q = similar(E)
    δp = similar(E)
    δq = similar(E)
    R = similar(E)
    Z = similar(E)
    λ = 1000
    iterations = 1
    w = 0.25*[0 1 0;1 0 1;0 1 0];
    x,y = 1:N,1:M
    wx = (2 .* π .* x) ./ M;
    wy = (2 .* π .* y) ./ N;
    for i = 1:iterations
        percent = (i/iterations)*100
        @show i,percent
        p_ = imfilter(p,reflect(w),"replicate")
        q_ = imfilter(q,reflect(w),"replicate")
        R = (ρ.*(-I[1] .* p .- I[2].* q .+ I[3]))./sqrt.(1 .+ p.^2 .+ q.^2)
        pq = (1 .+ p.^2 .+ q.^2)
        dR_dp = ((-ρ * I[1]) ./ sqrt.(pq)) .+ ((-I[1] * ρ) .* p .- I[2] * ρ .* q .+ I[3] * ρ) .* (-1 .* p .* (pq .^(-3/2)))
        dR_dq = ((-ρ*I[2]) ./ sqrt.(pq)) .+ (-I[1] * ρ .* p - I[2] * ρ .* q .+ I[3] * ρ) .* (-1 .* q .* (pq .^(-3/2)))
        p = p_ .+ (1/(4*λ))*(E.-R).*dR_dp;
        q = q_ .+ (1/(4*λ))*(E.-R).*dR_dq;
        Cp = fft(p)
        Cq = fft(q)
        C = -1im.*(wx .* Cp .+ wy .* Cq)./(wx.^2 .+ wy.^2)
        Z = abs.(ifft(C))
        p = ifft(1im .* wx .* C);
        q = ifft(1im .* wy .* C);
    end
    return Z
end
