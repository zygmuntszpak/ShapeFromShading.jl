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
            C[i] = (-1im * wx[i] * p[i] + 1im * wy[i] * q[i]) / (2*π*(wx[i]^2 + wy[i]^2 + eps()))
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

# function solve_EulerLagrange2(ρ, I, iterations, p, q, R, λ, E, Z)
#     bound = (2:(first(size(R))-1),(2:last(size(R))-1))
#     for i = 1:iterations
#         pᵏ⁺¹ = copy(p)
#         qᵏ⁺¹ = copy(q)
#         for i in CartesianIndices(bound)
#             #calculate reflectance map
#             R[i] = (ρ * (-I[1] * p[i] - I[2] * q[i] + I[3])) / sqrt(1 + p[i]^2 + q[i]^2)
#             pq = (1 + p[i]^2 + q[i]^2)
#             dRδp = ((-ρ * I[1]) / sqrt(pq)) + ((-I[1] * ρ) * p[i] - I[2] * ρ * q[i] + I[3] * ρ) * (-1 * p[i] * (pq^(-3/2)))
#             dRδq = ((-ρ * I[2]) / sqrt(pq)) + (-I[1] * ρ * p[i] - I[2] * ρ * q[i] + I[3] * ρ) * (-1 * q[i] * (pq^(-3/2)))
#             pbar = (p[i[1]+1, i[2]] + p[i[1]-1, i[2]]) / 2
#             qbar = (q[i[1], i[2]+1] + q[i[1], i[2]-1]) / 2
#             pt = (p[i[1]+1, i[2]+1] + p[i[1]-1, i[2]-1] - p[i[1]+1, i[2]-1] - p[i[1]-1, i[2]+1]) / 4
#             qt = (q[i[1]+1, i[2]+1] + q[i[1]-1, i[2]-1] - q[i[1]+1, i[2]-1] - q[i[1]-1, i[2]+1]) / 4
#             pᵏ⁺¹[i] = pbar - qt/2 + (1 / (2 * λ)) * (E[i] - R[i]) * dRδq
#             qᵏ⁺¹[i] = qbar - pt/2 + (1 / (2 * λ)) * (E[i] - R[i]) * dRδp
#         end
#         p = pᵏ⁺¹
#         q = qᵏ⁺¹
#     end
#     return p, q
# end


# WIP only works for odd dims
function setup_transform_values(N, M)
    wx = [j for i in 1:M, j in -0.5:(1/(N-1)):0.5]
    wy = [j for j in -0.5:(1/(M-1)):0.5, i in 1:N]
    wx = ifftshift(wx)
    wy = ifftshift(wy)
    return wx, wy
end

function setup_transform_values_pentland(M, N)
    wx = zeros(M, N)
    wy = zeros(M, N)
    for i in CartesianIndices(wx)
        wx[i] = (2 * π * i[2]) / M
        wy[i] = (2 * π * (N - i[1] + 1)) / N
    end
    return wx, wy
end

# Sets up xy grid as array
function setup_xy(dim)
    x = zeros(dim , dim)
    y = zeros(dim , dim)
    xyrange = ceil(Int, -dim/2):floor(Int, (dim/2))
    for i in CartesianIndices(x)
        x[i] = xyrange[i[2]]
        y[i] = -xyrange[i[1]]
    end
    return x,y
end

function gen_mask(p::AbstractArray, q::AbstractArray,σ::Real)
    pq, pp = imgradients(p, KernelFactors.prewitt)
    qq, qp = imgradients(q, KernelFactors.prewitt)
    # Violation of Schwarz integrability condition + large jump in one derivative without a corrisponding jump in the other
    mask = abs.(pq - qp) + abs.(pp - qq)
    t = mean(mask)+σ*std(mask)
    for i in CartesianIndices(mask)
          mask[i] = mask[i] <= t ? 1 : 0
    end
    stop = false
    while !stop
        stop = true
        for i in CartesianIndices(mask)
            if mask[i] != 0 && i[1] != 1 && i[2] != 1 && i[1] != first(size(p)) && i[2] != last(size(p))
                count = 0
                count += mask[i[1]+1,i[2]] == 0 ? 1 : 0
                count += mask[i[1]-1,i[2]] == 0 ? 1 : 0
                count += mask[i[1],i[2]+1] == 0 ? 1 : 0
                count += mask[i[1],i[2]-1] == 0 ? 1 : 0
                if count > 2
                    mask[i] = 0
                    stop = false
                end
            end
        end
    end
    splitMask = copy(mask)
    while findfirst(x -> x==1 , splitMask) != nothing
        # @show count(x -> x == 1, splitMask)
        seed = [findfirst(x -> x==1 , splitMask)]
        val = maximum(splitMask) + 1
        while length(seed) > 0
            i = 1
            # @show length(seed)
            if splitMask[seed[i]] == 1
                splitMask[seed[i]] = val
                if seed[i][1] != 1
                    if splitMask[seed[i][1]-1, seed[i][2]] == 1
                        push!(seed, CartesianIndex(seed[i][1]-1, seed[i][2]))
                    end
                end
                if seed[i][2] != last(size(mask))
                    if splitMask[seed[i][1], seed[i][2]+1] == 1
                        push!(seed, CartesianIndex(seed[i][1], seed[i][2]+1))
                    end
                end
                if seed[i][1] != first(size(mask))
                    if splitMask[seed[i][1]+1, seed[i][2]] == 1
                        push!(seed, CartesianIndex(seed[i][1]+1, seed[i][2]))
                    end
                end
                if seed[i][2] != 1
                    if splitMask[seed[i][1], seed[i][2]-1] == 1
                        push!(seed, CartesianIndex(seed[i][1], seed[i][2]-1))
                    end
                end
            end
            deleteat!(seed,1)
            # @show length(seed)
        end
    end
    layerMask = zeros(Float64, size(mask)..., Int(maximum(splitMask)))
    layerMask[:,:,1] = mask
    for i in CartesianIndices(splitMask)
        if splitMask[i] != 0
            layerMask[i,Int(splitMask[i])] = 1
        end
    end
    return layerMask
end
