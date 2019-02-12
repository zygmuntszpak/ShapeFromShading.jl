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
    return Z
end

function setup_transform_values(M, N)
    wx = zeros(M, N)
    wy = zeros(M, N)
    for i in CartesianIndices(wx)
        wx[i] = (2 * π * i[2]) / M
        wy[i] = (2 * π * i[1]) / N
    end
    return wx, wy
end
