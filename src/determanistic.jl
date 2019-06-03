@doc raw"""
WIP
"""
function retrieve_surface(algorithm::Determanistic, img::AbstractArray, p::AbstractArray, q::AbstractArray, integration_factor::Real = 10, smoothness::Real = 20)
    E = Array{Float64}(img)
    Eₘ = maximum(E)
    N, M = size(E)
    λᵢ = integration_factor
    λₛ = smoothness
    R = zeros(Float64,axes(E))
    Z = zeros(axes(E))
    Δϵ = zeros(Float64,(size(E)...,2))
    val = sqrt(2 * length(E))*2 + 10
    prevd = 0.0001
    while val > sqrt(2 * length(E))
        ϵ = 0.0
        for i = 1:(M-1)
            for j = 1:(N-1)
                ϵ += (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i])^2 + λᵢ*((p[j,i+1]-p[j,i])-(q[j+1,i]-q[j,i]))^2 + λₛ*((p[j,i+1]-p[j,i])^2 + (p[j+1,i]-p[j,i])^2 + (q[j,i+1]-q[j,i])^2 + (q[j+1,i]-q[j,i])^2)
                Δϵ[j,i,1] = -2 * Eₘ * p[j,i] * (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i]) * (1/(1+p[j,i]^2+q[j,i]^2)^(3/2)) - 2*λᵢ*((p[j+1,i]-p[j,i])-(q[j,i+1]-q[j,i])) - 2*λₛ*((p[j,i+1]-p[j,i]) + (p[j+1,i]-p[j,i]))
                Δϵ[j,i,2] = -2 * Eₘ * q[j,i] * (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i]) * (1/(1+p[j,i]^2+q[j,i]^2)^(3/2)) + 2*λᵢ*((p[j+1,i]-p[j,i])-(q[j,i+1]-q[j,i])) - 2*λₛ*((q[j,i+1]-q[j,i]) + (q[j+1,i]-q[j,i]))
            end
        end
        val = norm(Δϵ)
        dmin = 1
        ϵmin = typemax(Float64)
        for d = prevd/10:prevd/10:prevd*10
            ϵtemp = 0.0
            ptemp = copy(p)
            qtemp = copy(q)
            for i = 1:(M-1)
                for j = 1:(N-1)
                    ptemp[j,i] = p[j,i] - d*Δϵ[j,i,1]
                    qtemp[j,i] = q[j,i] - d*Δϵ[j,i,2]
                end
            end
            for i = 1:(M-1)
                for j = 1:(N-1)
                    ϵtemp += (Eₘ/sqrt(1+(ptemp[j,i])^2+qtemp[j,i]^2)-E[j,i])^2 + λᵢ*((ptemp[j,i+1]-ptemp[j,i])-(qtemp[j+1,i]-qtemp[j,i]))^2 + λₛ*((ptemp[j,i+1]-ptemp[j,i])^2 + (ptemp[j+1,i]-ptemp[j,i])^2 + (qtemp[j,i+1]-qtemp[j,i])^2 + (qtemp[j+1,i]-qtemp[j,i])^2)
                end
            end
            if ϵtemp < ϵmin
                ϵmin = ϵtemp
                dmin = d
            end
        end
        for i = 1:(M-1)
            for j = 1:(N-1)
                p[j,i] = p[j,i] - dmin*Δϵ[j,i,1]
                q[j,i] = q[j,i] - dmin*Δϵ[j,i,2]
            end
        end
        prevd = dmin
        @show maximum(Δϵ)
        @show norm(Δϵ), ϵ, dmin
    end
    return Z,p,q
end

function f(E, pq, λᵢ, λₛ)
    Eₘ = maximum(E)
    val = 0.0
    N, M = size(E)
    for i = 1:(M-1)
        for j = 1:(N-1)
            val += (Eₘ/sqrt(1+pq[j,i,1]^2+pq[j,i,2]^2)-E[j,i])^2 + λᵢ*((pq[j,i+1,1]-pq[j,i,1])-(pq[j+1,i,2]-pq[j,i,2]))^2 + λₛ*((pq[j,i+1,1]-pq[j,i,1])^2 + (pq[j+1,i,1]-pq[j,i,1])^2 + (pq[j,i+1,2]-pq[j,i,2])^2 + (pq[j+1,i,2]-pq[j,i,2])^2)
        end
    end
    return val
end

function retrieve_surface(algorithm::Determanistic2, img::AbstractArray, p::AbstractArray, q::AbstractArray, integration_factor::Real = 10, smoothness::Real = 20)

end
