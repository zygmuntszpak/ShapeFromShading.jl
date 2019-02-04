#solves the EulerLagrane equations setup in DiscreteShape and DiscreteShapeBound.
function solve_EulerLagrange(ρ, I, iterations, δp, δq, w, p, q, R, λ, wx, wy, E, Z)
    pq = zeros(Complex{Float64}, size(p))
    dRδp = zeros(Complex{Float64}, size(R))
    dRδq = zeros(Complex{Float64}, size(R))
    C = zeros(Complex{Float64}, size(R))
    for i = 1:iterations
        δp = imfilter(p, reflect(w), "replicate")
        δq = imfilter(q, reflect(w), "replicate")

        for i in CartesianIndices(R)
            R[i] = (ρ * (-I[1] * p[i] - I[2] * q[i] + I[3])) / sqrt(1 + p[i]^2 + q[i]^2)
            pq[i] = (1 + p[i]^2 + q[i]^2)
            dRδp[i] = ((-ρ * I[1]) / sqrt(pq[i])) + ((-I[1] * ρ) * p[i] - I[2] * ρ * q[i] + I[3] * ρ) * (-1 * p[i] * (pq[i]^(-3/2)))
            dRδq[i] = ((-ρ * I[2]) / sqrt(pq[i])) + (-I[1] * ρ * p[i] - I[2] * ρ * q[i] + I[3] * ρ) * (-1 * q[i] * (pq[i]^(-3/2)))
            p[i] = δp[i] + (1 / (4 * λ)) * (E[i] - R[i]) * dRδp[i]
            q[i] = δq[i] + (1 / (4 * λ)) * (E[i] - R[i]) * dRδq[i]
        end

        fft!(p)
        fft!(q)

        for i in CartesianIndices(C)
            C[i] = -1im * (wx[i] * p[i] + wy[i] * q[i]) / (wx[i]^2 + wy[i]^2)
            p[i] = 1im * wx[i] * C[i]
            q[i] = 1im * wy[i] * C[i]
        end

        ifft!(p)
        ifft!(q)
        ifft!(C)

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
