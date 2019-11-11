#converts gradient field into heightmap using Fourier method
function (scheme::Frankot)(pin::AbstractArray, qin::AbstractArray)
    p = Complex{Float64}.(pin)
    q = Complex{Float64}.(qin)
    wx, wy = setup_transform_values(last(size(p)),first(size(p)))
    fft!(p)
    fft!(q)
    Z = zeros(Complex{Float64}, size(p))
    Z = (-1im .* wx .* p .+ 1im .* wy .* q) ./ (2 .* π .* (wx.^2 .+ wy.^2 .+ eps()))
    ifft!(Z)
    Z = abs.(Z)
    return Z
end

#converts gradient field into heightmap using average of two integrals
function (scheme::Path)(p::AbstractArray, q::AbstractArray)
    R, C = size(p)
    Z₁ = zeros(Float64, size(p))
    Z₂ = zeros(Float64, size(p))

    #path 1
    for i = 2:R
        Z₁[i,1] = Z₁[i-1,1] - q[i,1]
    end
    for i = 1:R
        for j = 2:C
            Z₁[i,j] = Z₁[i,j-1] + p[i,j]
        end
    end

    #path 2
    for i = 2:C
        Z₂[1,i] = Z₂[1,i-1] + p[1,i]
    end
    for i = 2:R
        for j = 1:C
            Z₂[i,j] = Z₂[i-1,j] - q[i,j]
        end
    end
    Z = (Z₁ .+ Z₂) / 2
    return Z
end

#converts gradient field into heightmap using integral culculated from average gradient at point
function (scheme::SplitPath)(p::AbstractArray, q::AbstractArray)
    R, C = size(p)
    Z = zeros(Float64, size(p))
    for i = 2:R
        Z[i,1]=Z[i-1,1] + q[i,1]
        Z[1,i]=Z[1,i-1] - p[1,i]
    end
    for i = 2:R
        for j = i:R
            Z[i,j] = ((Z[i,j-1] + p[i,j]) + (Z[i-1,j] - q[i,j])) / 2
            Z[j,i] = ((Z[j,i-1] + p[j,i]) + (Z[j-1,i] - q[j,i])) / 2
        end
    end
    return Z
end

#Horn and Brooks
function (scheme::Horn)(p::AbstractArray, q::AbstractArray)#, iter::Real = 10000, ϵ::Real = 1)
    iter = copy(scheme.max_iter)
    ϵ = copy(scheme.ϵ)
    Z = zeros(Float64, size(p))
    Zᵏ⁺¹ = zeros(Float64, size(p))
    h = zeros(Float64, size(p))
    v = zeros(Float64, size(p))
    hv = zeros(Float64, size(p))
    R, C = size(p)
    for i = 1:R
        Z[i,1] = p[i,1]
        Z[i,R] = p[i,R]
    end
    for i = 1:C
        Z[1,i] = q[1,i]
        Z[C,i] = q[C,i]
    end
    for i = 2:(R-1)
        for j = 2:(C-1)
            h[i,j] = (p[i,j+1] - p[i,j-1]) / 2
            v[i,j] = (q[i+1,j] - q[i-1,j]) / 2
            hv[i,j] = (ϵ / 4) * (h[i,j] - v[i,j])
        end
    end
    for k = 1:iter
        copyto!(Zᵏ⁺¹, Z)
        for i = 2:(R-1)
            for j = 2:(C-1)
                Zᵏ⁺¹[i,j] = (Z[i-1,j] + Z[i+1,j] + Z[i,j-1] + Z[i,j+1]) / 4
                Zᵏ⁺¹[i,j] = Zᵏ⁺¹[i,j] - hv[i,j]
            end
        end
        Z = Zᵏ⁺¹
    end
    return Z
end

#Durou and Courteille
function (scheme::Durou)(p::AbstractArray, q::AbstractArray, iter::Real = 1000, ϵ::Real = 1)
    iter = copy(scheme.max_iter)
    ϵ = copy(scheme.ϵ)
    Z = zeros(Float64, size(p))
    Zᵏ⁺¹ = zeros(Float64, size(p))
    h = zeros(Float64, size(p))
    v = zeros(Float64, size(p))
    hv = zeros(Float64, size(p))
    R, C = size(p)
    for i = 1:(R-1)
        for j = 1:(C-1)
            h[i,j] = (p[i,j+1] + p[i,j]) / 2
            v[i,j] = (q[i+1,j] + q[i,j]) / 2
            hv[i,j] = (ϵ / 2) * (h[i,j] - v[i,j])
        end
    end
    for k = 1:iter
        copyto!(Zᵏ⁺¹, Z)
        for i = 1:(R-1)
            for j = 1:(C-1)
                Zᵏ⁺¹[i,j] = (Z[i+1,j] + Z[i,j+1]) / 2
                Zᵏ⁺¹[i,j] = Zᵏ⁺¹[i,j] - hv[i,j]
            end
        end
        Z = Zᵏ⁺¹
    end
    return Z
end

# Generates matricies Dᵤ⁺, Dᵤ⁻, Dᵥ⁺ and Dᵥ⁻ as per eq 11, 12, 13, 14 and 17
function gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻,T,mask)
    U,V = size(mask)
    for i in CartesianIndices(mask)
        #u+
        if i[1] + 1 <= U && mask[i] == 1 && mask[i[1]+1,i[2]] == 1
            Dᵤ⁺[i[1]+(i[2]-1)*U, i[1]+(i[2]-1)*U] = -1.0
            Dᵤ⁺[i[1]+(i[2]-1)*U, i[1]+(i[2]-1)*U+1] = 1.0
        end
        #u-
        if i[1] - 1 > 0 && mask[i] == 1 && mask[i[1]-1,i[2]] == 1
            Dᵤ⁻[i[1]+(i[2]-1)*U, i[1]+(i[2]-1)*U] = 1.0
            Dᵤ⁻[i[1]+(i[2]-1)*U, i[1]+(i[2]-1)*U-1] = -1.0
        end
        #v+
        if i[2] + 1 <= V && mask[i] == 1 && mask[i[1],i[2]+1] == 1
            Dᵥ⁺[i[1]+(i[2]-1)*V, i[1]+(i[2]-1)*V] = -1.0
            Dᵥ⁺[i[1]+(i[2]-1)*V, i[1]+(i[2])*V] = 1.0
        end
        #v-
        if i[2] - 1 > 0 && mask[i] == 1 && mask[i[1],i[2]-1] == 1
            Dᵥ⁻[i[1]+(i[2]-1)*V, i[1]+(i[2]-1)*V] = 1.0
            Dᵥ⁻[i[1]+(i[2]-1)*V, i[1]+(i[2]-2)*V] = -1.0
        end
    end
end

function (scheme::Quadratic)(pIn::AbstractArray, qIn::AbstractArray) #, λ::AbstractArray, z⁰::AbstractArray, maskIn::AbstractArray = fill(1.0, size(z⁰)))
    z⁰ = copy(scheme.z)
    λ = copy(scheme.λ)
    mask = copy(scheme.mask)
    z = copy(z⁰)
    mask = rotr90(mask)
    λ = rotr90(λ)
    p = rotr90(copy(pIn)).*mask
    q = rotr90(copy(qIn)).*mask
    T = first(size(z))
    index = CartesianIndices(z)
    index = reshape(index, length(index),1)

    Dᵤ⁺ = spzeros(Float64, T^2, T^2)
    Dᵤ⁻ = spzeros(Float64, T^2, T^2)
    Dᵥ⁺ = spzeros(Float64, T^2, T^2)
    Dᵥ⁻ = spzeros(Float64, T^2, T^2)
    gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻, T, mask)

    # Calculate Dᵤ, Dᵥ and L as per eq 23 and 24
    Dᵤ = 0.5*(transpose(Dᵤ⁺) .+ transpose(Dᵤ⁻))
    Dᵥ = 0.5*(transpose(Dᵥ⁺) .+ transpose(Dᵥ⁻))
    L = 0.5*((transpose(Dᵤ⁺) * Dᵤ⁺) .+ (transpose(Dᵤ⁻) * Dᵤ⁻) .+ (transpose(Dᵥ⁺) * Dᵥ⁺) .+ (transpose(Dᵥ⁻) * Dᵥ⁻))

    Λ = spzeros(Float64,T^2,T^2)
    A = spzeros(Float64,T^2,T^2)
    P = zeros(Float64,T^2)
    Q = zeros(Float64,T^2)
    Z⁰ = zeros(Float64,T^2)
    Z = zeros(Float64,T^2)

    # Vectorize inputs
    for j in CartesianIndices(z)
        Λ[j[1]+(j[2]-1)*T,j[1]+(j[2]-1)*T] = sqrt(λ[j])
        P[j[1]+(j[2]-1)*T] = p[j]
        Q[j[1]+(j[2]-1)*T] = q[j]
        Z⁰[j[1]+(j[2]-1)*T] = z⁰[j]
        Z[j[1]+(j[2]-1)*T] = z[j]
    end

    # Calculate A and b as per eq 23 abd 24
    A = L .+ (Λ^2)
    b = Dᵤ*P .+ Dᵥ*Q .+ (Λ^2)*Z⁰

    # Solve system
    pl = AMGPreconditioner{RugeStuben}(A)
    cg!(Z, A, b, tol=10.0^-4, Pl = pl)

    for j in eachindex(index)
        z[index[j]] = Z[j]
    end
    z = z.*mask
    z = rotl90(z)
    return z
end

function (scheme::TotalVariation)(pIn::AbstractArray, qIn::AbstractArray)
    z⁰ = copy(scheme.z)
    λ = copy(scheme.λ)
    mask = copy(scheme.mask)
    α = copy(scheme.α)
    max_iter = copy(scheme.max_iter)
    z = copy(z⁰)
    mask = rotr90(mask)
    λ = rotr90(λ)
    p = rotr90(copy(pIn)).*mask
    q = rotr90(copy(qIn)).*mask
    T = first(size(z))
    index = CartesianIndices(z)
    index = reshape(index, length(index),1)

    Dᵤ⁺ = spzeros(Float64,T^2,T^2)
    Dᵤ⁻ = spzeros(Float64,T^2,T^2)
    Dᵥ⁺ = spzeros(Float64,T^2,T^2)
    Dᵥ⁻ = spzeros(Float64,T^2,T^2)
    gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻, T, mask)
    Λ = spzeros(Float64,T^2,T^2)
    A = spzeros(Float64,T^2,T^2)
    P = zeros(Float64,T^2)
    Q = zeros(Float64,T^2)
    Z⁰ = zeros(Float64,T^2)
    Z = zeros(Float64,T^2)

    # Vectorize inputs
    for j in CartesianIndices(z)
        Λ[j[1]+(j[2]-1)*T,j[1]+(j[2]-1)*T] = sqrt(λ[j])
        P[j[1]+(j[2]-1)*T] = p[j]
        Q[j[1]+(j[2]-1)*T] = q[j]
        Z⁰[j[1]+(j[2]-1)*T] = z⁰[j]
        Z[j[1]+(j[2]-1)*T] = z[j]
    end

    # Calculate A as per eq 59
    A = (α/8)*(((transpose(Dᵤ⁺) * Dᵤ⁺) .+ (transpose(Dᵥ⁺) * Dᵥ⁺)) .+ ((transpose(Dᵤ⁻) * Dᵤ⁻) .+ (transpose(Dᵥ⁻) * Dᵥ⁻)) .+ ((transpose(Dᵤ⁻) * Dᵤ⁻) .+ (transpose(Dᵥ⁺) * Dᵥ⁺)) .+ ((transpose(Dᵤ⁺) * Dᵤ⁺) .+ (transpose(Dᵥ⁻) * Dᵥ⁻))) .+ (Λ^2)
    pl = AMGPreconditioner{RugeStuben}(A)

    # Calculate b, s, r, P and Q in each combination of {+,-} as per equations 57, 61, 62
    bᵤ⁺⁺ = zeros(Float64,T^2)
    bᵥ⁺⁺ = zeros(Float64,T^2)
    sᵤ⁺⁺ = Dᵤ⁺*Z - P + bᵤ⁺⁺
    sᵥ⁺⁺ = Dᵥ⁺*Z - Q + bᵥ⁺⁺
    rᵤ⁺⁺ = max(norm(sᵤ⁺⁺) - (4/α),0)*(sᵤ⁺⁺./(norm(sᵤ⁺⁺)+eps(Float64)))
    rᵥ⁺⁺ = max(norm(sᵥ⁺⁺) - (4/α),0)*(sᵥ⁺⁺./(norm(sᵥ⁺⁺)+eps(Float64)))
    bᵤ⁺⁺ = bᵤ⁺⁺ + Dᵤ⁺*Z - P - rᵤ⁺⁺
    bᵥ⁺⁺ = bᵥ⁺⁺ + Dᵥ⁺*Z - Q - rᵥ⁺⁺
    P⁺⁺ = P + rᵤ⁺⁺ - bᵤ⁺⁺
    Q⁺⁺ = Q + rᵥ⁺⁺ - bᵥ⁺⁺

    bᵤ⁺⁻ = zeros(Float64,T^2)
    bᵥ⁺⁻ = zeros(Float64,T^2)
    sᵤ⁺⁻ = Dᵤ⁺*Z - P + bᵤ⁺⁻
    sᵥ⁺⁻ = Dᵥ⁻*Z - Q + bᵥ⁺⁻
    rᵤ⁺⁻ = max(norm(sᵤ⁺⁻) - (4/α),0)*(sᵤ⁺⁻./(norm(sᵤ⁺⁻)+eps(Float64)))
    rᵥ⁺⁻ = max(norm(sᵥ⁺⁻) - (4/α),0)*(sᵥ⁺⁻./(norm(sᵥ⁺⁻)+eps(Float64)))
    bᵤ⁺⁻ = bᵤ⁺⁻ + Dᵤ⁺*Z - P - rᵤ⁺⁻
    bᵥ⁺⁻ = bᵥ⁺⁻ + Dᵥ⁻*Z - Q - rᵥ⁺⁻
    P⁺⁻ = P + rᵤ⁺⁻ - bᵤ⁺⁻
    Q⁺⁻ = Q + rᵥ⁺⁻ - bᵥ⁺⁻

    bᵤ⁻⁺ = zeros(Float64,T^2)
    bᵥ⁻⁺ = zeros(Float64,T^2)
    sᵤ⁻⁺ = Dᵤ⁻*Z - P + bᵤ⁻⁺
    sᵥ⁻⁺ = Dᵥ⁺*Z - Q + bᵥ⁻⁺
    rᵤ⁻⁺ = max(norm(sᵤ⁻⁺) - (4/α),0)*(sᵤ⁻⁺./(norm(sᵤ⁻⁺)+eps(Float64)))
    rᵥ⁻⁺ = max(norm(sᵥ⁻⁺) - (4/α),0)*(sᵥ⁻⁺./(norm(sᵥ⁻⁺)+eps(Float64)))
    bᵤ⁻⁺ = bᵤ⁻⁺ + Dᵤ⁻*Z - P - rᵤ⁻⁺
    bᵥ⁻⁺ = bᵥ⁻⁺ + Dᵥ⁺*Z - Q - rᵥ⁻⁺
    P⁻⁺ = P + rᵤ⁻⁺ - bᵤ⁻⁺
    Q⁻⁺ = Q + rᵥ⁻⁺ - bᵥ⁻⁺

    bᵤ⁻⁻ = zeros(Float64,T^2)
    bᵥ⁻⁻ = zeros(Float64,T^2)
    sᵤ⁻⁻ = Dᵤ⁻*Z - P + bᵤ⁻⁻
    sᵥ⁻⁻ = Dᵥ⁻*Z - Q + bᵥ⁻⁻
    rᵤ⁻⁻ = max(norm(sᵤ⁻⁻) - (4/α),0)*(sᵤ⁻⁻./(norm(sᵤ⁻⁻)+eps(Float64)))
    rᵥ⁻⁻ = max(norm(sᵥ⁻⁻) - (4/α),0)*(sᵥ⁻⁻./(norm(sᵥ⁻⁻)+eps(Float64)))
    bᵤ⁻⁻ = bᵤ⁻⁻ + Dᵤ⁻*Z - P - rᵤ⁻⁻
    bᵥ⁻⁻ = bᵥ⁻⁻ + Dᵥ⁻*Z - Q - rᵥ⁻⁻
    P⁻⁻ = P + rᵤ⁻⁻ - bᵤ⁻⁻
    Q⁻⁻ = Q + rᵥ⁻⁻ - bᵥ⁻⁻

    # Calculate bₜᵥ as per eq 60
    bₜᵥ = (transpose(Dᵤ⁺) * P⁺⁺) .+ (transpose(Dᵥ⁺) * Q⁺⁺) .+ (transpose(Dᵤ⁺) * P⁺⁻) .+ (transpose(Dᵥ⁻) * Q⁺⁻) .+ (transpose(Dᵤ⁻) * P⁻⁺) .+ (transpose(Dᵥ⁺) * Q⁻⁺) .+ (transpose(Dᵤ⁻) * P⁻⁻) .+ (transpose(Dᵥ⁻) * Q⁻⁻) .+ (Λ^2)*Z⁰
    for i = 1:max_iter
        cg!(Z, A, bₜᵥ, tol=10.0^-4, Pl = pl)

        # Recalculate b, s, r, P and Q in each combination of {+,-} as per equations 57, 61, 62
        sᵤ⁺⁺ = Dᵤ⁺*Z - P + bᵤ⁺⁺
        sᵥ⁺⁺ = Dᵥ⁺*Z - Q + bᵥ⁺⁺
        rᵤ⁺⁺ = max(norm(sᵤ⁺⁺) - (4/α),0)*(sᵤ⁺⁺./(norm(sᵤ⁺⁺)+eps(Float64)))
        rᵥ⁺⁺ = max(norm(sᵥ⁺⁺) - (4/α),0)*(sᵥ⁺⁺./(norm(sᵥ⁺⁺)+eps(Float64)))
        bᵤ⁺⁺ = bᵤ⁺⁺ + Dᵤ⁺*Z - P - rᵤ⁺⁺
        bᵥ⁺⁺ = bᵥ⁺⁺ + Dᵥ⁺*Z - Q - rᵥ⁺⁺
        P⁺⁺ = P + rᵤ⁺⁺ - bᵤ⁺⁺
        Q⁺⁺ = Q + rᵥ⁺⁺ - bᵥ⁺⁺

        sᵤ⁺⁻ = Dᵤ⁺*Z - P + bᵤ⁺⁻
        sᵥ⁺⁻ = Dᵥ⁻*Z - Q + bᵥ⁺⁻
        rᵤ⁺⁻ = max(norm(sᵤ⁺⁻) - (4/α),0)*(sᵤ⁺⁻./(norm(sᵤ⁺⁻)+eps(Float64)))
        rᵥ⁺⁻ = max(norm(sᵥ⁺⁻) - (4/α),0)*(sᵥ⁺⁻./(norm(sᵥ⁺⁻)+eps(Float64)))
        bᵤ⁺⁻ = bᵤ⁺⁻ + Dᵤ⁺*Z - P - rᵤ⁺⁻
        bᵥ⁺⁻ = bᵥ⁺⁻ + Dᵥ⁻*Z - Q - rᵥ⁺⁻
        P⁺⁻ = P + rᵤ⁺⁻ - bᵤ⁺⁻
        Q⁺⁻ = Q + rᵥ⁺⁻ - bᵥ⁺⁻

        sᵤ⁻⁺ = Dᵤ⁻*Z - P + bᵤ⁻⁺
        sᵥ⁻⁺ = Dᵥ⁺*Z - Q + bᵥ⁻⁺
        rᵤ⁻⁺ = max(norm(sᵤ⁻⁺) - (4/α),0)*(sᵤ⁻⁺./(norm(sᵤ⁻⁺)+eps(Float64)))
        rᵥ⁻⁺ = max(norm(sᵥ⁻⁺) - (4/α),0)*(sᵥ⁻⁺./(norm(sᵥ⁻⁺)+eps(Float64)))
        bᵤ⁻⁺ = bᵤ⁻⁺ + Dᵤ⁻*Z - P - rᵤ⁻⁺
        bᵥ⁻⁺ = bᵥ⁻⁺ + Dᵥ⁺*Z - Q - rᵥ⁻⁺
        P⁻⁺ = P + rᵤ⁻⁺ - bᵤ⁻⁺
        Q⁻⁺ = Q + rᵥ⁻⁺ - bᵥ⁻⁺

        sᵤ⁻⁻ = Dᵤ⁻*Z - P + bᵤ⁻⁻
        sᵥ⁻⁻ = Dᵥ⁻*Z - Q + bᵥ⁻⁻
        rᵤ⁻⁻ = max(norm(sᵤ⁻⁻) - (4/α),0)*(sᵤ⁻⁻./(norm(sᵤ⁻⁻)+eps(Float64)))
        rᵥ⁻⁻ = max(norm(sᵥ⁻⁻) - (4/α),0)*(sᵥ⁻⁻./(norm(sᵥ⁻⁻)+eps(Float64)))
        bᵤ⁻⁻ = bᵤ⁻⁻ + Dᵤ⁻*Z - P - rᵤ⁻⁻
        bᵥ⁻⁻ = bᵥ⁻⁻ + Dᵥ⁻*Z - Q - rᵥ⁻⁻
        P⁻⁻ = P + rᵤ⁻⁻ - bᵤ⁻⁻
        Q⁻⁻ = Q + rᵥ⁻⁻ - bᵥ⁻⁻

        # Recalculate bₜᵥ as per eq 60
        bₜᵥ = (transpose(Dᵤ⁺) * P⁺⁺) .+ (transpose(Dᵥ⁺) * Q⁺⁺) .+ (transpose(Dᵤ⁺) * P⁺⁻) .+ (transpose(Dᵥ⁻) * Q⁺⁻) .+ (transpose(Dᵤ⁻) * P⁻⁺) .+ (transpose(Dᵥ⁺) * Q⁻⁺) .+ (transpose(Dᵤ⁻) * P⁻⁻) .+ (transpose(Dᵥ⁻) * Q⁻⁻) .+ (Λ^2)*Z⁰
    end
    for j in eachindex(index)
        z[index[j]] = Z[j]
    end
    z = z.*mask
    z = rotl90(z)
    return z
end

function (scheme::NonConvex1)(pIn::AbstractArray, qIn::AbstractArray)
    z⁰ = copy(scheme.z)
    λ = copy(scheme.λ)
    mask = copy(scheme.mask)
    β = copy(scheme.β)
    max_iter = copy(scheme.max_iter)
    z = copy(z⁰)
    mask = rotr90(mask)
    λ = rotr90(λ)
    p = rotr90(copy(pIn)).*mask
    q = rotr90(copy(qIn)).*mask
    T = first(size(z))
    index = CartesianIndices(z)
    index = reshape(index, length(index),1)

    Dᵤ⁺ = spzeros(Float64,T^2,T^2)
    Dᵤ⁻ = spzeros(Float64,T^2,T^2)
    Dᵥ⁺ = spzeros(Float64,T^2,T^2)
    Dᵥ⁻ = spzeros(Float64,T^2,T^2)
    gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻, T, mask)
    Λ = zeros(Float64,T^2,T^2)
    P = zeros(Float64,T^2)
    Q = zeros(Float64,T^2)
    Z⁰ = zeros(Float64,T^2)
    Z = zeros(Float64,T^2)

    # Vectorize inputs
    for j in CartesianIndices(z)
        Λ[j[1]+(j[2]-1)*T,j[1]+(j[2]-1)*T] = sqrt(λ[j])
        P[j[1]+(j[2]-1)*T] = p[j]
        Q[j[1]+(j[2]-1)*T] = q[j]
        Z⁰[j[1]+(j[2]-1)*T] = z⁰[j]
        Z[j[1]+(j[2]-1)*T] = z[j]
    end
    Λ = Diagonal(Λ)
    Zold = copy(Z)
    Ztemp = copy(Z)
    α₁ = 0.8
    α₂ = 0.8

    # Initialize Lₙ (Lipschitz constant) such that lazy backtracking starts with α₁ = α₂.
    # This is mostly arbitary and value will approach needed value over time.
    Lₙ = (18.0/α₁)*(1-α₂)
    η = 1.1
    i = 1
    while i <= max_iter && (norm(Z - Zold) > 10.0^-4.0 || i == 1)
        # Calculate f₁, ∇f₁ and previous f₁ (f₁ⁿ⁻¹ )
        ∇f₁ = ((transpose(Dᵤ⁺)*(Dᵤ⁺*Z-P))/(norm(Dᵤ⁺*Z-P)^2+β^2)) + ((transpose(Dᵥ⁺)*(Dᵥ⁺*Z-Q))/(norm(Dᵥ⁺*Z-Q)^2+β^2)) + ((transpose(Dᵤ⁻)*(Dᵤ⁻*Z-P))/(norm(Dᵤ⁻*Z-P)^2+β^2)) + ((transpose(Dᵥ⁻)*(Dᵥ⁻*Z-Q))/(norm(Dᵥ⁻*Z-Q)^2+β^2))
        f₁ⁿ⁻¹ = log(norm(Dᵤ⁺*Zold-P)^2+β^2) + log(norm(Dᵥ⁺*Zold-Q)^2+β^2) + log(norm(Dᵤ⁻*Zold-P)^2+β^2) + log(norm(Dᵥ⁻*Zold-Q)^2+β^2)
        f₁ = log(norm(Dᵤ⁺*Z-P)^2+β^2) + log(norm(Dᵥ⁺*Z-Q)^2+β^2) + log(norm(Dᵤ⁻*Z-P)^2+β^2) + log(norm(Dᵥ⁻*Z-Q)^2+β^2)

        # Calculate lazy backtracking condition as per ipiano algortihm 4
        while f₁ > f₁ⁿ⁻¹  + dot(∇f₁, Z - Zold) + (Lₙ/2)*norm(Z - Zold)
            Lₙ = Lₙ*η
        end

        # Calculate new α₁ as per ipiano algortihm 4
        α₁ = (2*(1-α₂)/Lₙ)*0.9

        # Calcualte x̂ as per eq 69 then calculate new Z as per eq 71
        x̂ = Z - α₁.*∇f₁ + α₂.*(Z-Zold)
        Ztemp = inv(I+(2*α₁).*Λ^2)*(x̂+(2*α₁).*Λ*Z⁰)
        copyto!(Zold, Z)
        copyto!(Z, Ztemp)
        i += 1
    end
    println("Converged after ", i-1, " iterations")
    for j in eachindex(index)
        z[index[j]] = Z[j]
    end
    z = rotl90(z)
    return z
end

function (scheme::AnisotropicDiffusion)(pIn::AbstractArray, qIn::AbstractArray) #, λ::AbstractArray, z⁰::AbstractArray, maskIn::AbstractArray = fill(1.0, size(z⁰)), μ::Real = 1.0, ν::Real = 10.0; max_iter::Int = 1000)
    z⁰ = copy(scheme.z)
    λ = copy(scheme.λ)
    mask = copy(scheme.mask)
    ν = copy(scheme.ν)
    μ = copy(scheme.μ)
    max_iter = copy(scheme.max_iter)
    z = copy(z⁰)
    mask = rotr90(mask)
    λ = rotr90(λ)
    p = rotr90(copy(pIn)).*mask
    q = rotr90(copy(qIn)).*mask
    T = first(size(z))
    index = CartesianIndices(z)
    index = reshape(index, length(index),1)

    Dᵤ⁺ = spzeros(Float64,T^2,T^2)
    Dᵤ⁻ = spzeros(Float64,T^2,T^2)
    Dᵥ⁺ = spzeros(Float64,T^2,T^2)
    Dᵥ⁻ = spzeros(Float64,T^2,T^2)
    gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻,T,mask)
    Λ = zeros(Float64,T^2,T^2)
    P = zeros(Float64,T^2)
    Q = zeros(Float64,T^2)
    Z⁰ = zeros(Float64,T^2)
    Z = zeros(Float64,T^2)

    # Vectorize inputs
    for j in CartesianIndices(z)
        Λ[j[1]+(j[2]-1)*T,j[1]+(j[2]-1)*T] = sqrt(λ[j])
        P[j[1]+(j[2]-1)*T] = p[j]
        Q[j[1]+(j[2]-1)*T] = q[j]
        Z⁰[j[1]+(j[2]-1)*T] = z⁰[j]
        Z[j[1]+(j[2]-1)*T] = z[j]
    end
    Λ = sparse(Diagonal(Λ))

    # Calcualte A and B for each conination in {+,-} as per eq 90 and 91
    A⁺⁺ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
    A⁻⁺ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
    A⁺⁻ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))
    A⁻⁻ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))
    B⁺⁺ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
    B⁻⁺ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
    B⁺⁻ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))
    B⁻⁻ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))

    # Calculate A and b as derived from eq 92 in form Az=b
    # A = A⁺⁺*transpose(A⁺⁺)*Dᵤ⁺*transpose(Dᵤ⁺) + B⁺⁺*transpose(B⁺⁺)*Dᵥ⁺*transpose(Dᵥ⁺) + A⁺⁻*transpose(A⁺⁻)*Dᵤ⁺*transpose(Dᵤ⁺) + B⁺⁻*transpose(B⁺⁻)*Dᵥ⁻*transpose(Dᵥ⁻) + A⁻⁺*transpose(A⁻⁺)*Dᵤ⁻*transpose(Dᵤ⁻) + B⁻⁺*transpose(B⁻⁺)*Dᵥ⁺*transpose(Dᵥ⁺) + A⁻⁻*transpose(A⁻⁻)*Dᵤ⁻*transpose(Dᵤ⁻) + B⁻⁻*transpose(B⁻⁻)*Dᵥ⁻*transpose(Dᵥ⁻)
    A = (A⁺⁺*Dᵤ⁺)*transpose(A⁺⁺*Dᵤ⁺) + (B⁺⁺*Dᵥ⁺)*transpose(B⁺⁺*Dᵥ⁺) + (A⁺⁻*Dᵤ⁺)*transpose(A⁺⁻*Dᵤ⁺) + (B⁺⁻*Dᵥ⁻)*transpose(B⁺⁻*Dᵥ⁻) + (A⁻⁺*Dᵤ⁻)*transpose( A⁻⁺*Dᵤ⁻) + (B⁻⁺*Dᵥ⁺)*transpose(B⁻⁺*Dᵥ⁺) + (A⁻⁻*Dᵤ⁻)*transpose(A⁻⁻*Dᵤ⁻) + (B⁻⁻*Dᵥ⁻)*transpose(B⁻⁻*Dᵥ⁻)
    A = A/4.0 + Λ^2
    # b = A⁺⁺*transpose(A⁺⁺)*transpose(Dᵤ⁺)*P + B⁺⁺*transpose(B⁺⁺)*transpose(Dᵥ⁺)*Q + A⁺⁻*transpose(A⁺⁻)*transpose(Dᵤ⁺)*P + B⁺⁻*transpose(B⁺⁻)*transpose(Dᵥ⁻)*Q + A⁻⁺*transpose(A⁻⁺)*transpose(Dᵤ⁻)*P + B⁻⁺*transpose(B⁻⁺)*transpose(Dᵥ⁺)*Q + A⁻⁻*transpose(A⁻⁻)*transpose(Dᵤ⁻)*P + B⁻⁻*transpose(B⁻⁻)*transpose(Dᵥ⁻)*Q
    b = (transpose(A⁺⁺*Dᵤ⁺)*A⁺⁺)*P + (transpose(B⁺⁺*Dᵥ⁺)*B⁺⁺)*Q + (transpose(A⁺⁻*Dᵤ⁺)*A⁺⁻)*P + (transpose(B⁺⁻*Dᵥ⁻)*B⁺⁻)*Q + (transpose(A⁻⁺*Dᵤ⁻)*A⁻⁺)*P + (transpose(B⁻⁺*Dᵥ⁺)*B⁻⁺)*Q + (transpose(A⁻⁻*Dᵤ⁻)*A⁻⁻)*P + (transpose(B⁻⁻*Dᵥ⁻)*B⁻⁻)*Q
    b = b/4.0 + (Λ^2)*Z⁰
    A₀ = factorize(A)
    for i = 1:max_iter
        # Solve system at iteration k
        Z = A₀\b

        # Recalcualte A and B for each conination in {+,-} as per eq 90 and 91
        A⁺⁺ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
        A⁻⁺ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
        A⁺⁻ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))
        A⁻⁻ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))
        B⁺⁺ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
        B⁻⁺ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
        B⁺⁻ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))
        B⁻⁻ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))

        # Recalculate A and b as derived from eq 92 in form Az=b
        # A = A⁺⁺*transpose(A⁺⁺)*Dᵤ⁺*transpose(Dᵤ⁺) + B⁺⁺*transpose(B⁺⁺)*Dᵥ⁺*transpose(Dᵥ⁺) + A⁺⁻*transpose(A⁺⁻)*Dᵤ⁺*transpose(Dᵤ⁺) + B⁺⁻*transpose(B⁺⁻)*Dᵥ⁻*transpose(Dᵥ⁻) + A⁻⁺*transpose(A⁻⁺)*Dᵤ⁻*transpose(Dᵤ⁻) + B⁻⁺*transpose(B⁻⁺)*Dᵥ⁺*transpose(Dᵥ⁺) + A⁻⁻*transpose(A⁻⁻)*Dᵤ⁻*transpose(Dᵤ⁻) + B⁻⁻*transpose(B⁻⁻)*Dᵥ⁻*transpose(Dᵥ⁻)
        A = (A⁺⁺*Dᵤ⁺)*transpose(A⁺⁺*Dᵤ⁺) + (B⁺⁺*Dᵥ⁺)*transpose(B⁺⁺*Dᵥ⁺) + (A⁺⁻*Dᵤ⁺)*transpose(A⁺⁻*Dᵤ⁺) + (B⁺⁻*Dᵥ⁻)*transpose(B⁺⁻*Dᵥ⁻) + (A⁻⁺*Dᵤ⁻)*transpose( A⁻⁺*Dᵤ⁻) + (B⁻⁺*Dᵥ⁺)*transpose(B⁻⁺*Dᵥ⁺) + (A⁻⁻*Dᵤ⁻)*transpose(A⁻⁻*Dᵤ⁻) + (B⁻⁻*Dᵥ⁻)*transpose(B⁻⁻*Dᵥ⁻)
        A = A/4.0 + Λ^2
        # b = A⁺⁺*transpose(A⁺⁺)*transpose(Dᵤ⁺)*P + B⁺⁺*transpose(B⁺⁺)*transpose(Dᵥ⁺)*Q + A⁺⁻*transpose(A⁺⁻)*transpose(Dᵤ⁺)*P + B⁺⁻*transpose(B⁺⁻)*transpose(Dᵥ⁻)*Q + A⁻⁺*transpose(A⁻⁺)*transpose(Dᵤ⁻)*P + B⁻⁺*transpose(B⁻⁺)*transpose(Dᵥ⁺)*Q + A⁻⁻*transpose(A⁻⁻)*transpose(Dᵤ⁻)*P + B⁻⁻*transpose(B⁻⁻)*transpose(Dᵥ⁻)*Q
        b = (transpose(A⁺⁺*Dᵤ⁺)*A⁺⁺)*P + (transpose(B⁺⁺*Dᵥ⁺)*B⁺⁺)*Q + (transpose(A⁺⁻*Dᵤ⁺)*A⁺⁻)*P + (transpose(B⁺⁻*Dᵥ⁻)*B⁺⁻)*Q + (transpose(A⁻⁺*Dᵤ⁻)*A⁻⁺)*P + (transpose(B⁻⁺*Dᵥ⁺)*B⁻⁺)*Q + (transpose(A⁻⁻*Dᵤ⁻)*A⁻⁻)*P + (transpose(B⁻⁻*Dᵥ⁻)*B⁻⁻)*Q
        b = b/4.0 + (Λ^2)*Z⁰
        A₀ = factorize(A)
    end

    for j in eachindex(index)
        z[index[j]] = Z[j]
    end
    z = z.*mask
    z = rotl90(z)
    return z
end

function (scheme::MumfordShah)(pIn::AbstractArray, qIn::AbstractArray)#, λ::AbstractArray, z⁰::AbstractArray, maskIn::AbstractArray = fill(1.0, size(z⁰)), μ::Real = 1.0, ϵ::Real = 1.0; max_iter::Int = 1000)
    z⁰ = copy(scheme.z)
    # z⁰ = copy(z⁰)
    # mask = copy(maskIn)
    λ = copy(scheme.λ)
    mask = copy(scheme.mask)
    ϵ = copy(scheme.ϵ)
    μ = copy(scheme.μ)
    max_iter = copy(scheme.max_iter)
    z = copy(z⁰)
    mask = rotr90(mask)
    λ = rotr90(copy(λ))
    p = rotr90(copy(pIn)).*mask
    q = rotr90(copy(qIn)).*mask
    T = first(size(z))
    index = CartesianIndices(z)
    index = reshape(index, length(index),1)

    Dᵤ⁺ = spzeros(Float64,T^2,T^2)
    Dᵤ⁻ = spzeros(Float64,T^2,T^2)
    Dᵥ⁺ = spzeros(Float64,T^2,T^2)
    Dᵥ⁻ = spzeros(Float64,T^2,T^2)
    gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻,T,mask)
    Λ = zeros(Float64,T^2,T^2)
    P = zeros(Float64,T^2)
    Q = zeros(Float64,T^2)
    Z⁰ = zeros(Float64,T^2)
    Z = zeros(Float64,T^2)
    # Vectorize inputs
    for j in CartesianIndices(z)
        Λ[j[1]+(j[2]-1)*T,j[1]+(j[2]-1)*T] = sqrt(λ[j])
        P[j[1]+(j[2]-1)*T] = p[j]
        Q[j[1]+(j[2]-1)*T] = q[j]
        Z⁰[j[1]+(j[2]-1)*T] = z⁰[j]
        Z[j[1]+(j[2]-1)*T] = z[j]
    end
    Λ = sparse(Diagonal(Λ))
    wᵤ⁺ = fill(1.0, T^2)
    wᵤ⁻ = fill(1.0, T^2)
    wᵥ⁺ = fill(1.0, T^2)
    wᵥ⁻ = fill(1.0, T^2)

    Wᵤ⁺ = Diagonal(wᵤ⁺)
    Wᵤ⁻ = Diagonal(wᵤ⁻)
    Wᵥ⁺ = Diagonal(wᵥ⁺)
    Wᵥ⁻ = Diagonal(wᵥ⁻)

    A1 = μ*(transpose(Dᵤ⁺)*transpose(Wᵤ⁺)*Wᵤ⁺*Dᵤ⁺ + transpose(Dᵤ⁻)*transpose(Wᵤ⁻)*Wᵤ⁻*Dᵤ⁻ + transpose(Dᵥ⁺)*transpose(Wᵥ⁺)*Wᵥ⁺*Dᵥ⁺ + transpose(Dᵥ⁻)*transpose(Wᵥ⁻)*Wᵥ⁻*Dᵥ⁻)
    A1 = A1 + Λ^2
    b = μ*((transpose(Dᵤ⁺)*transpose(Wᵤ⁺)*Wᵤ⁺)*P + (transpose(Dᵤ⁻)*transpose(Wᵤ⁻)*Wᵤ⁻)*P + (transpose(Dᵥ⁺)*transpose(Wᵥ⁺)*Wᵥ⁺)*Q + (transpose(Dᵥ⁻)*transpose(Wᵥ⁻)*Wᵥ⁻)*Q)
    b = b + (Λ^2)*Z⁰

    b₀ =  Diagonal(fill(1.0/(4.0*ϵ), T^2))
    DDᵤ⁺ = (ϵ*transpose(Dᵤ⁺)*Dᵤ⁺+b₀)
    DDᵤ⁻ = (ϵ*transpose(Dᵤ⁻)*Dᵤ⁻+b₀)
    DDᵥ⁺ = (ϵ*transpose(Dᵥ⁺)*Dᵥ⁺+b₀)
    DDᵥ⁻ = (ϵ*transpose(Dᵥ⁻)*Dᵥ⁻+b₀)
    A2 = ((Diagonal(Dᵤ⁻*Z-P)^2)*μ) + DDᵤ⁻
    b2 = fill(1.0/(4.0*ϵ), T^2)
    for i = 1:max_iter
        cg!(Z, A1, b)

        A2 = ((Diagonal(Dᵤ⁺*Z-P)^2)*μ) + DDᵤ⁺
        cg!(wᵤ⁺, A2, b2)
        Wᵤ⁺ = Diagonal(wᵤ⁺)

        A2 = ((Diagonal(Dᵤ⁻*Z-P)^2)*μ) + DDᵤ⁻
        cg!(wᵤ⁻, A2, b2)
        Wᵤ⁻ = Diagonal(wᵤ⁻)

        A2 = ((Diagonal(Dᵥ⁺*Z-Q)^2)*μ) + DDᵥ⁺
        cg!(wᵥ⁺, A2, b2)
        Wᵥ⁺ = Diagonal(wᵥ⁺)

        A2 = ((Diagonal(Dᵥ⁻*Z-Q)^2)*μ) + DDᵥ⁻
        cg!(wᵥ⁻, A2, b2)
        Wᵥ⁻ = Diagonal(wᵥ⁻)

        A1 = μ*(transpose(Dᵤ⁺)*transpose(Wᵤ⁺)*Wᵤ⁺*Dᵤ⁺ + transpose(Dᵤ⁻)*transpose(Wᵤ⁻)*Wᵤ⁻*Dᵤ⁻ + transpose(Dᵥ⁺)*transpose(Wᵥ⁺)*Wᵥ⁺*Dᵥ⁺ + transpose(Dᵥ⁻)*transpose(Wᵥ⁻)*Wᵥ⁻*Dᵥ⁻)
        A1 = A1 + Λ^2
        b = μ*((transpose(Dᵤ⁺)*transpose(Wᵤ⁺)*Wᵤ⁺)*P + (transpose(Dᵤ⁻)*transpose(Wᵤ⁻)*Wᵤ⁻)*P + (transpose(Dᵥ⁺)*transpose(Wᵥ⁺)*Wᵥ⁺)*Q + (transpose(Dᵥ⁻)*transpose(Wᵥ⁻)*Wᵥ⁻)*Q)
        b = b + (Λ^2)*Z⁰
    end

    for j in eachindex(index)
        z[index[j]] = Z[j]
    end
    z = z.*mask
    z = rotl90(z)
    return z
end
