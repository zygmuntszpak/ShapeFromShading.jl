#solves the EulerLagrane equations setup in DiscreteShape and DiscreteShapeBound.
function solve_EulerLagrange(ρ, I, iterations, p, q, R, λ, E, Z)
    #initilize common variables
    M, N = size(E)
    δp = zeros(Complex{Float64},axes(E))
    δq = zeros(Complex{Float64},axes(E))
    pq = zeros(Complex{Float64}, size(p))
    dRδp = zeros(Complex{Float64}, size(R))
    dRδq = zeros(Complex{Float64}, size(R))
    C = zeros(Complex{Float64}, size(R))
    w = centered(0.25*[0.0 1.0 0.0;1.0 0.0 1.0;0.0 1.0 0.0])
    wx, wy = setup_transform_values(M, N)
    for i = 1:iterations
        #calculate second derivatives of p and q
        δp = imfilter(p, reflect(w), "replicate")
        δq = imfilter(q, reflect(w), "replicate")

        for i in CartesianIndices(R)
            #calculate reflectance map
            R[i] = (ρ * (-I[1] * p[i] - I[2] * q[i] + I[3])) / sqrt(1 + p[i]^2
                + q[i]^2)
            #calculate partial derivatives of p and q in respect to p and q
            pq[i] = (1 + p[i]^2 + q[i]^2)
            dRδp[i] = ((-ρ * I[1]) / sqrt(pq[i])) + ((-I[1] * ρ) * p[i] - I[2]
                * ρ * q[i] + I[3] * ρ) * (-1 * p[i] * (pq[i]^(-3/2)))

            dRδq[i] = ((-ρ * I[2]) / sqrt(pq[i])) + (-I[1] * ρ * p[i] - I[2]
                * ρ * q[i] + I[3] * ρ) * (-1 * q[i] * (pq[i]^(-3/2)))
            #calculate new surface normals
            p[i] = δp[i] + (1 / (4 * λ)) * (E[i] - R[i]) * dRδp[i]
            q[i] = δq[i] + (1 / (4 * λ)) * (E[i] - R[i]) * dRδq[i]
        end

        #calculate Transform of surface normals
        fft!(p)
        fft!(q)

        for i in CartesianIndices(C)
            #compute C and prepare to update p and q
            C[i] = -1im * (wx[i] * p[i] + wy[i] * q[i]) / (wx[i]^2 + wy[i]^2)
            p[i] = 1im * wx[i] * C[i]
            q[i] = 1im * wy[i] * C[i]
        end

        #update p and q
        ifft!(p)
        ifft!(q)

        #recover Z by taking inverse transform of C
        ifft!(C)

        #take abs of C to get Z
        for i in CartesianIndices(Z)
            Z[i] = abs(C[i])
        end
    end
    return Z, p, q
end

function setup_transform_values(M, N)
    wx = zeros(M, N)
    wy = zeros(M, N)
    for i in CartesianIndices(wx)
        wx[i] = (2 * π * i[2]) / M
        wy[i] = (2 * π * (N - i[1] + 1)) / N
    end
    return wx, wy
end

#converts gradient field into heightmap using Fourier method
function convert_gradient(pIn::AbstractArray, qIn::AbstractArray)
    p = Complex{Float64}.(pIn)
    q = Complex{Float64}.(qIn)
    wx, wy = setup_transform_values(last(size(p)),first(size(p)))
    Cp = fft(p)
    Cq = fft(q)
    Z = zeros(Complex{Float64}, size(p))
    for i in CartesianIndices(Z)
        Z[i] = -1im * (wx[i] * Cp[i] + wy[i] * Cq[i]) / (wx[i]^2 + wy[i]^2)
    end
    ifft!(Z)
    Z = abs.(Z)
    return Z
end

#converts gradient field into heightmap using average of two integrals
function convert_gradient2(p::AbstractArray, q::AbstractArray)
    R, C = size(p)
    Z₁ = zeros(Float64, size(p))
    Z₂ = zeros(Float64, size(p))

    #path 1
    for i = 2:R
        Z₁[i,1] = Z₁[i-1,1] + q[i,1]
    end
    for i = 1:R
        for j = 2:C
            Z₁[i,j] = Z₁[i,j-1] - p[i,j]
        end
    end

    #path 2
    for i = 2:C
        Z₂[1,i] = Z₂[1,i-1] - p[1,i]
    end
    for i = 2:R
        for j = 1:C
            Z₂[i,j] = Z₂[i-1,j] + q[i,j]
        end
    end
    Z = (Z₁ .+ Z₂) / 2
    return Z
end

#converts gradient field into heightmap using integral culculated from average gradient at point
function convert_gradient3(p::AbstractArray, q::AbstractArray)
    R, C = size(p)
    Z = zeros(Float64, size(p))
    for i = 2:R
        Z[i,1]=Z[i-1,1] + q[i,1]
        Z[1,i]=Z[1,i-1] - p[1,i]
    end
    for i = 2:R
        for j = i:R
            Z[i,j] = ((Z[i,j-1] - p[i,j]) + (Z[i-1,j] + q[i,j])) / 2
            Z[j,i] = ((Z[j,i-1] - p[j,i]) + (Z[j-1,i] + q[j,i])) / 2
        end
    end
    return Z
end

#converts gradient field into heightmap using weighted average of number of integrals from each corner
function convert_gradient4(p::AbstractArray, q::AbstractArray)
    R, C = size(p)
    Z = zeros(Float64, R, C)
    Z₁ = zeros(Float64, R, C)
    Z₂ = zeros(Float64, R, C)
    Z₃ = zeros(Float64, R, C)
    Z₄ = zeros(Float64, R, C)
    Z₅ = zeros(Float64, R, C)
    Z₆ = zeros(Float64, R, C)
    Z₇ = zeros(Float64, R, C)
    Z₈ = zeros(Float64, R, C)

    #top left path 1
    for i = 2:R
        Z₁[i,1] = Z₁[i-1,1] + q[i,1]
    end
    for i = 1:R
        for j = 2:C
            Z₁[i,j] = Z₁[i,j-1] - p[i,j]
        end
    end

    #path 2
    for i = 2:C
        Z₂[1,i] = Z₂[1,i-1] - p[1,i]
    end
    for i = 2:R
        for j = 1:C
            Z₂[i,j] = Z₂[i-1,j] + q[i,j]
        end
    end

    #bottem left path 3
    for i = R-1:-1:1
        Z₃[i,1] = Z₃[i+1,1] - q[i,1]
    end
    for i = R:-1:1
        for j = 2:C
            Z₃[i,j] = Z₃[i,j-1] - p[i,j]
        end
    end

    #path 4
    for i = 2:C
        Z₄[R,i] = Z₄[R,i-1] - p[R,i]
    end
    for i = R-1:-1:1
        for j = 1:C
            Z₄[i,j] = Z₄[i+1,j] - q[i,j]
        end
    end

    #top Right path 5
    for i = 2:R
        Z₅[i,C] = Z₅[i-1,C] + q[i,C]
    end
    for i = 1:R
        for j = C-1:-1:1
            Z₅[i,j] = Z₅[i,j+1] + p[i,j]
        end
    end

    #path 6
    for i = C-1:-1:1
        Z₆[1,i] = Z₆[1,i+1] + p[1,i]
    end
    for i = 2:R
        for j = C:-1:1
            Z₆[i,j] = Z₆[i-1,j] + q[i,j]
        end
    end

    #bottem right path 7
    for i = R-1:-1:1
        Z₇[i,C] = Z₇[i+1,C] - q[i,C]
    end
    for i = R:-1:1
        for j = C-1:-1:1
            Z₇[i,j] = Z₇[i,j+1] + p[i,j]
        end
    end

    #path 8
    for i = C-1:-1:1
        Z₈[R,i] = Z₈[R,i+1] + p[R,i]
    end
    for i = R-1:-1:1
        for j = C:-1:1
            Z₈[i,j] = Z₈[i+1,j] - q[i,j]
        end
    end

    #average
    for i = 1:R
        for j = 1:C
            dTL = sqrt(i^2 + j^2)
            dBL = sqrt((R - i + 1)^2 + j^2)
            dTR = sqrt(i^2 + (C - j + 1)^2)
            dBR = sqrt((R - i + 1)^2 + (C - j + 1)^2)
            ZTL = (Z₁[i,j] + Z₂[i,j]) / 2
            ZBL = (Z₃[i,j] + Z₄[i,j]) / 2
            ZTR = (Z₅[i,j] + Z₆[i,j]) / 2
            ZBR = (Z₇[i,j] + Z₈[i,j]) / 2
            Z[i,j] = (ZTL / dTL) + (ZBL / dBL) + (ZTR / dTR) + (ZBR / dBR)
            Z[i,j] = Z[i,j] / 4
        end
    end
    return Z
end
