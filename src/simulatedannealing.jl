function F(E, p, q, λᵢ, λₛ, Eₑ, Eₛ, Eᵢ, energy, add, j, i)
    Eₘ = maximum(E)
    ϵ = 0.0
    N, M = size(E)
    if add == true
        push!(Eₑ, 0.0)
        push!(Eₛ, 0.0)
        push!(Eᵢ, 0.0)
    end

    for k = (i-1):(i+1)
        for r = (j-1):(j+1)
            if k != N && r != M && k != 1 && r != 1
                # @show r,k
                # Optimize this!!!
                s = CartesianIndex(r,k)
                s1 = CartesianIndex(r+1,k)
                s2 = CartesianIndex(r,k+1)
                s3 = CartesianIndex(r-1,k)
                s4 = CartesianIndex(r,k-1)
                s5 = CartesianIndex(r+1,k-1)
                s6 = CartesianIndex(r-1,k+1)
                ϵ₀ = (Eₘ/sqrt(1+p[s]^2+q[s]^2)-E[s])^2
                ϵₛ = λₛ*((p[s1] - p[s])^2 + (q[s1] - q[s])^2 + (p[s2] - p[s])^2 + (q[s2] - q[s])^2 + (p[s3] - p[s])^2 + (q[s3] - q[s])^2 + (p[s4] - p[s])^2 + (q[s4] - q[s])^2)
                ϵᵢ = λᵢ*(((p[s1] - p[s]) - (q[s2] - q[s]))^2 + ((p[s] - p[s3]) - (q[s6] - q[s3]))^2 + ((p[s5] - p[s4]) - (q[s] - q[s4]))^2)
                ϵ += ϵ₀ + ϵᵢ + ϵₛ
                # @show ϵ₀, ϵᵢ, ϵₛ
                energy[i,j] = ϵ₀ + ϵᵢ + ϵₛ
                # energy[i,j] = ϵ₀
                if add == true
                    Eₑ[end] += ϵ₀
                    Eₛ[end] += ϵₛ
                    Eᵢ[end] += ϵᵢ
                end
                # energy[i,j] = ϵ₀
                # energy[i,j] = ϵᵢ
                # energy[i,j] = ϵₛ
            end
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

function (algorithm::SimulatedAnnealing)(img::AbstractArray)

    pin = algorithm.pin
    qin = algorithm.qin
    T₀ = algorithm.T₀
    α = algorithm.α
    integration_factor = algorithm.integration_factor
    smoothness = algorithm.smoothness
    max_iter = algorithm.max_iter
    debug = algorithm.debug

    ϵₑ = zeros(Float64, 0)
    ϵᵢ = zeros(Float64, 0)
    ϵₛ = zeros(Float64, 0)

    # Setup initial variables
    E = Float64.(copy(img)).*255
    Eₘ = maximum(E)
    p = copy(pin)
    q = copy(qin)
    pmax = maximum(abs.(p))
    qmax = maximum(abs.(q))
    ρₘ = sqrt(maximum(p.^2) + maximum(q.^2))
    pₙ = copy(pin)
    qₙ = copy(qin)
    p₀ = copy(pin)
    q₀ = copy(qin)
    λᵢ = integration_factor
    λₛ = smoothness
    k = 1
    Tₖ = T₀
    energy = zeros(size(E))
    ϵ = F(E, p, q, λᵢ, λₛ, ϵₑ, ϵₛ, ϵᵢ, energy, true, 20, 20)
    ϵ₀ = F(E, p, q, λᵢ, λₛ, ϵₑ, ϵₛ, ϵᵢ, energy, false, 20, 20)
    N, M = size(E)
    r2 = 1.0:1:32
    if debug == true
        Z = convert_gradient(Horn(),p,q)
        img = generate_surface(p, q, 1, [0,0,1], radius = 50, img_size=32, noise = 0, background = false)
        display(AbstractPlotting.PlotDisplay(), Makie.vbox(surface(r2, r2, Z), image(img'), heatmap(energy), arrows(r2, r2, pin[end:-1:1,:]', qin[end:-1:1,:]')))
    end

    while k < max_iter
        for i in CartesianIndices(E)
            if i[1] != N && i[2] != M && i[1] != 1 && i[2] != 1
                ϵ = F(E, p, q, λᵢ, λₛ, ϵₑ, ϵᵢ, ϵₛ, energy, false, i[1], i[2])
                θ = rand(0:10^-6:π/2)
                Φ = rand(0:10^-6:2*π)
                # ρₙ = tan(θ)
                ρₙ = sqrt(((Eₘ^2)/(E[i]^2)) - 1)
                if ρₙ > 2*ρₘ
                    ρₙ = 2*ρₘ
                end
                pₙ[i] = ρₙ*cos(Φ)
                qₙ[i] = ρₙ*sin(Φ)
                ϵₙ = F(E, pₙ, qₙ, λᵢ, λₛ, ϵₑ, ϵₛ, ϵᵢ, energy, false, i[1], i[2])
                ρᵢ = sqrt(p[i]^2 + q[i]^2)
                # ρₙ = sqrt(pₙ[i]^2 + qₙ[i]^2)
                ρ = (ρₙ * (1 + ρᵢ^2)^(3/2)) / (ρᵢ * (1 + ρₙ^2)^(3/2))
                R = P(ϵ, ϵₙ, ρ, Tₖ)
                prob = rand(0:0.00001:1)
                if prob < R
                    copyto!(p,pₙ)
                    copyto!(q,qₙ)
                else
                    copyto!(pₙ,p)
                    copyto!(qₙ,q)
                end
            end
        end
        ϵ = F(E, pₙ, qₙ, λᵢ, λₛ, ϵₑ, ϵₛ, ϵᵢ, energy, true, 20, 20)
        if k % 1000 == 0
            @show k
        end
        if k % 200 == 0 && debug == true
            r2 = 1.0:1:32
            Z = convert_gradient(Horn(),p,q)
            img = generate_surface(p, q, 1, [0,0,1], radius = 50, img_size=32, noise = 0, background = false)
            display(Makie.vbox(surface(r2, r2, Z), image(img'), heatmap((energy./maximum(abs.(energy)))), arrows(r2, r2, -p', q')))
            display(AbstractPlotting.PlotDisplay(), Makie.vbox(surface(r2, r2, Z), image(img'), heatmap((energy./maximum(abs.(energy)))), arrows(r2, r2, p[end:-1:1,:]', q[end:-1:1,:]')))
            # return energy, 0
        end
        if k % 100 == 0 && debug == true
            @show ϵ, Tₖ, k#, sum(energy)
        end
        k += 1
        Tₖ = (α^k)*T₀
    end
    if debug == true
        @show ϵ₀
    end
    return p,q,ϵₑ,ϵₛ,ϵᵢ
end
