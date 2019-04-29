#converts gradient field into heightmap using Fourier method
function convert_gradient(Scheme::Frankot, pIn::AbstractArray, qIn::AbstractArray)
    p = Complex{Float64}.(pIn)
    q = Complex{Float64}.(qIn)
    wx, wy = setup_transform_values(last(size(p)),first(size(p)))
    fft!(p)
    fft!(q)
    Z = zeros(Complex{Float64}, size(p))
    Z = (-1im .* wx .* p .+ 1im .* wy .* q) ./ (2 .* π .* (wx.^2 .+ wy.^2 .+ eps()))
    ifft!(Z)
    Z = abs.(Z)
    return Z
end

# function convert_gradient(Scheme::Frankot, pIn::AbstractArray, qIn::AbstractArray)
#     p = Complex{Float64}.(pIn)
#     q = Complex{Float64}.(qIn)
#     wx, wy = setup_transform_values(last(size(p)),first(size(p)))
#     Cp = fft(p)
#     Cq = fft(q)
#     Z = zeros(Complex{Float64}, size(p))
#     for i in CartesianIndices(Z)
#         Z[i] = -1im * (wx[i] * Cp[i] + wy[i] * Cq[i]) / (2*π*(wx[i]^2 + wy[i]^2))
#     end
#     ifft!(Z)
#     Z = abs.(Z)
#     return Z
# end

#converts gradient field into heightmap using average of two integrals
function convert_gradient(Scheme::Path, p::AbstractArray, q::AbstractArray)
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
function convert_gradient(Scheme::SplitPath, p::AbstractArray, q::AbstractArray)
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
function convert_gradient(Scheme::Horn, p::AbstractArray, q::AbstractArray, iter::Real = 10000, ϵ::Real = 1)
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
function convert_gradient(Scheme::Durou, p::AbstractArray, q::AbstractArray, iter::Real = 1000, ϵ::Real = 1)
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
