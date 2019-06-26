@doc raw"""
WIP
"""
function retrieve_surface(algorithm::Union{DeterministicCustom,Deterministic}, starting_shape::AbstractSyntheticShape, img::AbstractArray, integration_factor::Real = 10, smoothness::Real = 20; noise_val::Real = 0, radius=50)
    # Generate initial gradient and dipatch on desired function.
    p,q = synthetic_gradient(starting_shape, radius = radius, noise = noise_val)
    p,q = retrieve_surface(algorithm, img, p, q, integration_factor, smoothness)
    return p,q
end

function retrieve_surface(algorithm::DeterministicCustom, img::AbstractArray, pin::AbstractArray, qin::AbstractArray, integration_factor::Real = 10, smoothness::Real = 20)
    # Setup initial variables
    E = Array{Float64}(img)
    Eₘ = maximum(E)
    p = copy(pin)
    q = copy(qin)
    N, M = size(E)
    λᵢ = integration_factor
    λₛ = smoothness
    R = zeros(Float64,axes(E))
    Z = zeros(axes(E))
    Δϵ = zeros(Float64,(size(E)...,2))
    val = sqrt(2 * length(E))*2 + 10
    prevd = 0.0001

    # Loop untill norm of Δϵ is less then desired value
    while val > sqrt(2 * length(E))

        # Calculate the energy
        ϵ = 0.0
        for i = 1:(M-1)
            for j = 1:(N-1)
                ϵ += (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i])^2 + λᵢ*((p[j,i+1]-p[j,i])-(q[j+1,i]-q[j,i]))^2 + λₛ*((p[j,i+1]-p[j,i])^2 + (p[j+1,i]-p[j,i])^2 + (q[j,i+1]-q[j,i])^2 + (q[j+1,i]-q[j,i])^2)
                if i>1 && j>1
                    Δϵ[j,i,1] = -2 * Eₘ * p[j,i] * (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i]) * (1/(1+p[j,i]^2+q[j,i]^2)^(3/2)) - 2*λᵢ*((p[j+1,i]-p[j,i])-(q[j,i+1]-q[j,i])) - 2*λₛ*((p[j,i+1]-p[j,i]) + (p[j+1,i]-p[j,i])) + 2*λᵢ*((p[j,i]-p[j-1,i])-(q[j-1,i+1]-q[j-1,i])) + 2*λₛ*((p[j,i]-p[j,i-1]) + (p[j,i]-p[j-1,i]))
                    Δϵ[j,i,2] = -2 * Eₘ * q[j,i] * (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i]) * (1/(1+p[j,i]^2+q[j,i]^2)^(3/2)) + 2*λᵢ*((p[j+1,i]-p[j,i])-(q[j,i+1]-q[j,i])) - 2*λₛ*((q[j,i+1]-q[j,i]) + (q[j+1,i]-q[j,i])) - 2*λᵢ*((p[j+1,i-1]-p[j,i-1])-(q[j,i]-q[j,i-1])) + 2*λₛ*((q[j,i]-q[j,i-1]) + (q[j,i]-q[j-1,i]))
                elseif i>1
                    Δϵ[j,i,1] = -2 * Eₘ * p[j,i] * (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i]) * (1/(1+p[j,i]^2+q[j,i]^2)^(3/2)) - 2*λᵢ*((p[j+1,i]-p[j,i])-(q[j,i+1]-q[j,i])) - 2*λₛ*((p[j,i+1]-p[j,i]) + (p[j+1,i]-p[j,i])) + 2*λₛ*(p[j,i]-p[j,i-1])
                    Δϵ[j,i,2] = -2 * Eₘ * q[j,i] * (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i]) * (1/(1+p[j,i]^2+q[j,i]^2)^(3/2)) + 2*λᵢ*((p[j+1,i]-p[j,i])-(q[j,i+1]-q[j,i])) - 2*λₛ*((q[j,i+1]-q[j,i]) + (q[j+1,i]-q[j,i])) - 2*λᵢ*((p[j+1,i-1]-p[j,i-1])-(q[j,i]-q[j,i-1])) + 2*λₛ*(q[j,i]-q[j,i-1])
                elseif j>1
                    Δϵ[j,i,1] = -2 * Eₘ * p[j,i] * (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i]) * (1/(1+p[j,i]^2+q[j,i]^2)^(3/2)) - 2*λᵢ*((p[j+1,i]-p[j,i])-(q[j,i+1]-q[j,i])) - 2*λₛ*((p[j,i+1]-p[j,i]) + (p[j+1,i]-p[j,i])) + 2*λᵢ*((p[j,i]-p[j-1,i])-(q[j-1,i+1]-q[j-1,i])) +2*λₛ*(p[j,i]-p[j-1,i])
                    Δϵ[j,i,2] = -2 * Eₘ * q[j,i] * (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i]) * (1/(1+p[j,i]^2+q[j,i]^2)^(3/2)) + 2*λᵢ*((p[j+1,i]-p[j,i])-(q[j,i+1]-q[j,i])) - 2*λₛ*((q[j,i+1]-q[j,i]) + (q[j+1,i]-q[j,i])) + 2*λᵢ*(q[j,i]-q[j-1,i])
                else
                    Δϵ[j,i,1] = -2 * Eₘ * p[j,i] * (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i]) * (1/(1+p[j,i]^2+q[j,i]^2)^(3/2)) - 2*λᵢ*((p[j+1,i]-p[j,i])-(q[j,i+1]-q[j,i])) - 2*λₛ*((p[j,i+1]-p[j,i]) + (p[j+1,i]-p[j,i]))
                    Δϵ[j,i,2] = -2 * Eₘ * q[j,i] * (Eₘ/sqrt(1+p[j,i]^2+q[j,i]^2)-E[j,i]) * (1/(1+p[j,i]^2+q[j,i]^2)^(3/2)) + 2*λᵢ*((p[j+1,i]-p[j,i])-(q[j,i+1]-q[j,i])) - 2*λₛ*((q[j,i+1]-q[j,i]) + (q[j+1,i]-q[j,i]))
                end
            end
        end

        # Update norm value
        val = norm(Δϵ)

        # Line search to find minimizor of ϵ along the direction of gradient
        dmin = 1
        ϵmin = typemax(Float64)
        for d = prevd/10:prevd/10:prevd*10
            # Setup new gradient field
            ϵtemp = 0.0
            ptemp = copy(p)
            qtemp = copy(q)
            for i = 1:(M-1)
                for j = 1:(N-1)
                    ptemp[j,i] = p[j,i] - d*Δϵ[j,i,1]
                    qtemp[j,i] = q[j,i] - d*Δϵ[j,i,2]
                end
            end

            # Calculate new energy
            for i = 1:(M-1)
                for j = 1:(N-1)
                    ϵtemp += (Eₘ/sqrt(1+(ptemp[j,i])^2+qtemp[j,i]^2)-E[j,i])^2 + λᵢ*((ptemp[j,i+1]-ptemp[j,i])-(qtemp[j+1,i]-qtemp[j,i]))^2 + λₛ*((ptemp[j,i+1]-ptemp[j,i])^2 + (ptemp[j+1,i]-ptemp[j,i])^2 + (qtemp[j,i+1]-qtemp[j,i])^2 + (qtemp[j+1,i]-qtemp[j,i])^2)
                end
            end

            # Check if new energy is less then the minimum found
            if ϵtemp < ϵmin
                ϵmin = ϵtemp
                dmin = d
            end
        end

        # Set gradients to minimize energy
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
    return p,q
end

# Energy function
function f(E, pq, λᵢ, λₛ)
    Eₘ = maximum(E)
    val = 0.0
    N, M = size(E)
    for i = 1:(M-1)
        for j = 1:(N-1)
            val += (Eₘ/sqrt(1+pq[j,i,1]^2+pq[j,i,2]^2)-E[j,i])^2 + λᵢ*((pq[j,i+1,1]-pq[j,i,1])-(pq[j+1,i,2]-pq[j,i,1]))^2 + λₛ*((pq[j,i+1,1]-pq[j,i,1])^2 + (pq[j+1,i,1]-pq[j,i,1])^2 + (pq[j,i+1,2]-pq[j,i,2])^2 + (pq[j+1,i,2]-pq[j,i,2])^2)
        end
    end
    return val
end

function retrieve_surface(algorithm::Deterministic, img::AbstractArray, pin::AbstractArray, qin::AbstractArray, integration_factor::Real = 10, smoothness::Real = 20)
    # Setup initial variables
    E = Float64.(img)
    x₀ = zeros(Float64,first(size(img)),last(size(img)),2)
    x₀[:,:,1] .= pin
    x₀[:,:,2] .= qin

    # Calcualte derivatives
    od = OnceDifferentiable(pq -> f(E,pq,integration_factor,smoothness), x₀; autodiff = :forward);
    a = zeros(size(x₀))

    # Run gradient descent
    res = optimize(od, x₀, GradientDescent(), Optim.Options(show_trace = true, g_tol = sqrt(2*length(E))))

    # Get minimised gradients
    out = Optim.minimizer(res)
    p = out[:,:,1]
    q = out[:,:,2]
    return p,q
end
