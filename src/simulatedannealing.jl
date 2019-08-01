function F(E, p, q, λᵢ, λₛ)
    Eₘ = maximum(E)
    ϵ = 0.0
    N, M = size(E)
    for i = 1:(M-1)
        for j = 1:(N-1)
            ϵ += (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i])^2 + λᵢ*((p[j,i+1]-p[j,i])-(q[j+1,i]-q[j,i]))^2 + λₛ*((p[j,i+1]-p[j,i])^2 + (p[j+1,i]-p[j,i])^2 + (q[j,i+1]-q[j,i])^2 + (q[j+1,i]-q[j,i])^2)
        end
    end
    return ϵ
end

function P(ϵ, ϵₙ, ρ, T)
    R = 1.0
    if ϵ < ϵₙ
        R = exp(-(ϵₙ-ϵ)/T)
        # @show R
        R = R*(ρ^(1/T-1))
        # R = 0.0
    end
    return R
end

function retrieve_surface(algorithm::SimulatedAnnealing, img::AbstractArray, pin::AbstractArray, qin::AbstractArray, T₀::Real = 500.0, α::Real = 0.99, integration_factor::Real = 10, smoothness::Real = 20)
    # Setup initial variables
    E = Float64.(img)
    p = copy(pin)
    q = copy(qin)
    pmax = maximum(p)
    qmax = maximum(q)
    pmin = minimum(p)
    qmin = minimum(q)
    pₙ = copy(pin)
    qₙ = copy(qin)
    p₀ = copy(pin)
    q₀ = copy(qin)
    λᵢ = integration_factor
    λₛ = smoothness
    k = 1
    Tₖ = T₀
    ϵ = F(E, p, q, λᵢ, λₛ)
    ϵ₀ = F(E, p, q, λᵢ, λₛ)
    while Tₖ > 0.0001
        for i in CartesianIndices(E)
            if i != CartesianIndex(size(E))
                ϵ = F(E, p, q, λᵢ, λₛ)
                pₙ[i] += rand(Normal(0,1))
                qₙ[i] += rand(Normal(0,1))
                pₙ[i] = max(pmin, min(pmax, pₙ[i]))
                qₙ[i] = max(qmin, min(qmax, qₙ[i]))
                ϵₙ = F(E, pₙ, qₙ, λᵢ, λₛ)
                ρᵢ = sqrt(p[i]^2 + q[i]^2)
                ρₙ = sqrt(pₙ[i]^2 + qₙ[i]^2)
                ρ = (ρₙ * (1 + ρᵢ^2)^(3/2)) / (ρᵢ * (1 + ρₙ^2)^(3/2))
                R = P(ϵ, ϵₙ, ρ, Tₖ)
                # @show R, ϵ, ϵₙ
                prob = rand(0:0.00001:1)
                if prob < R
                    p = copy(pₙ)
                    q = copy(qₙ)
                    # @show R, prob
                else
                    pₙ = copy(p)
                    qₙ = copy(q)
                end
            end
        end
        # p[first(size(E)),last(size(p))] = 0.0
        # q[first(size(E)),last(size(p))] = 0.0
        # if ϵ < ϵ₀
        #     ϵ₀ = ϵ
        #     copyto!(p₀, p)
        #     copyto!(q₀, q)
        # end
        if k % 100 == 0
            r2 = 1.0:1:32
            Z = convert_gradient(Horn(),p,q)
            img = generate_surface(p, q, 1, [0,0,1], radius = 50, img_size=32, noise = 0, background = false)
            display(AbstractPlotting.PlotDisplay(), Makie.vbox(surface(r2, r2, Z), image(img)))
        end
        @show ϵ, Tₖ, k
        k += 1
        Tₖ = (α^k)*T₀
    end
    @show ϵ₀
    return p,q
end
