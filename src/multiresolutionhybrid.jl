#WIP currently runs simulated SimulatedAnnealing in compressed img
function (algorithm::MultiResolutionHybrid)(img::AbstractArray)
    @warn "This function is WIP and will only return zeros."
    # T₀ = algorithm.T₀
    # α = algorithm.α
    # integration_factor = algorithm.integration_factor
    # smoothness = algorithm.smoothness
    # integration_factor2 = algorithm.integration_factor2
    # smoothness2 = algorithm.smoothness2
    # img_size = first(size(img))
    # E = Float64.(img)
    # Eₘ = maximum(E)
    # ρ = sqrt.((Eₘ^2)./(E.^2) .- 1)
    # ρ = imfilter(ρ, Kernel.gaussian(3))
    # E₁ = Eₘ ./ (1 .+ (ρ[1:2:end, 1:2:end,].^2))
    # ρ = sqrt.((Eₘ^2)./(E₁.^2) .- 1)
    # ρ = imfilter(ρ, Kernel.gaussian(3))
    # E₂ = Eₘ ./ (1 .+ (ρ[1:2:end, 1:2:end,].^2))
    # ρ = sqrt.((Eₘ^2)./(E₂.^2) .- 1)
    # ρ = imfilter(ρ, Kernel.gaussian(3))
    # E₃ = Eₘ ./ (1 .+ (ρ[1:2:end, 1:2:end,].^2))
    # p = rand(Float64, size(E₃))
    # q = rand(Float64, size(E₃))
    # pn, qn = SimulatedAnnealing(pin = p, qin = q, T₀ = T₀, α = α, integration_factor = integration_factor, smoothness = smoothness, max_iter = 1)(E₃)
    p = zeros(Float64, size(img))
    q = zeros(Float64, size(img))
    return p, q
end
