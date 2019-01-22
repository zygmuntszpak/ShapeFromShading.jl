function retrieve_surface(algorithm::Shah, img::AbstractArray, iterations::Int=200)
    ρ,I,σ,τ = estimate_img_properties(img)
    E = Float64.(img)
    E = E.*255
    M,N=size(E)
    p = zeros(axes(E))
    q = zeros(axes(E))
    Z = zeros(axes(E))
    Z_x = zeros(axes(E))
    Z_y = zeros(axes(E))
    ix = cos(τ) * tan(σ)
    iy = sin(τ) * tan(σ)
    R = zeros(axes(E))
    @inbounds for i = 1:iterations
        R = (cos(σ) .+ p .* cos(τ)*sin(σ) .+ q .* sin(τ)*sin(σ))./sqrt.(1 .+ p.^2 + q.^2)
        R = max.(0,R)
        f = E .- R
        df_dZ =(p+q).*(ix.*p + iy.*q .+ 1)./(sqrt.((1 .+ p.^2 + q.^2).^3)*sqrt(1 + ix^2 + iy^2))-(ix+iy)./(sqrt.(1 .+ p.^2 + q.^2)*sqrt(1 + ix^2 + iy^2))
        Z = Z - f./(df_dZ .+ eps())
        Z_x[2:M,:] = Z[1:M-1,:]
        Z_y[:,2:N] = Z[:,1:N-1]
        p = Z - Z_x;
        q = Z - Z_y;
    end
    Z = mapwindow(median!, abs.(Z), (21,21))
    return Z
end
