@doc raw"""
```
DeterministicCustom(p::AbstractArray, q::AbstractArray, integration_factor::Real = 10, smoothness::Real = 20, β::Real = 1.0, integration_scheme::AbstractIntegrator = Horn())
```
Same as [`Deterministic`](@ref) exept impliments a custom gradient descent and
line search algorithm.
# Example
The following example demonstrates the use of the `DeterministicCustom` method.
```julia
using ShapeFromShading, Makie, Images

#generate synthetic image and initial solution
img = generate_surface(SynthSphere(38))
pin, qin = synthetic_gradient(SynthSphere(38))

#calculate the heightmap
deterministicCustom = DeterministicCustom(p=pin,q=qin, integration_factor = 10, smoothness = 50, β = 1.0, integration_scheme = Horn())
Z,p,q = deterministicCustom(img)

#normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(abs.(Z))

#display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
surface(r, r, Z)
```
# Reference
[1] A. Crouzil, X. Descombes and J. Durou, "A multiresolution approach for shape from shading coupling deterministic and stochastic optimization," in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 25, no. 11, pp. 1416-1421, Nov. 2003. [doi: 10.1109/TPAMI.2003.1240116](https://doi.org/10.1109/TPAMI.2003.1240116 )
"""
function (algorithm::DeterministicCustom)(img::AbstractArray)
    # Setup initial variables
    E = Array{Float64}(img)
    Eₘ = maximum(E)
    p = copy(algorithm.p)
    q = copy(algorithm.q)
    N, M = size(E)
    λᵢ = algorithm.integration_factor
    λₛ = algorithm.smoothness
    β = algorithm.β
    R = zeros(Float64,axes(E))
    Z = zeros(axes(E))
    Δϵ = zeros(Float64,(size(E)...,2))
    val = sqrt(2 * length(E))*2 + 10
    prevd = 0.0001

    # Loop untill norm of Δϵ is less then desired value
    while val > β*sqrt(2 * length(E))

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
    Z = algorithm.integration_scheme(p,q)
    return Z,p,q
end

# Energy function
function f(E, pq, λᵢ, λₛ)
    Eₘ = maximum(E)
    val = 0.0
    N, M = size(E)
    for j = 1:(M-1)
        for i = 1:(N-1)
            val += (Eₘ/sqrt(1+pq[j,i,1]^2+pq[j,i,2]^2)-E[j,i])^2 + λᵢ*((pq[j,i+1,1]-pq[j,i,1])-(pq[j+1,i,2]-pq[j,i,2]))^2 + λₛ*((pq[j,i+1,1]-pq[j,i,1])^2 + (pq[j+1,i,1]-pq[j,i,1])^2 + (pq[j,i+1,2]-pq[j,i,2])^2 + (pq[j+1,i,2]-pq[j,i,2])^2)
        end
    end
    return val
end

@doc raw"""
```
Deterministic(p::AbstractArray, q::AbstractArray, integration_factor::Real = 10, smoothness::Real = 20, β::Real = 1.0, integration_scheme::AbstractIntegrator = Horn())
```
Impliments a deterministic method for reconstructing the surface using a gradient
descent method. This version makes use of a [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl )
solver.
# Output
Returns a `Deterministic` functor which can be called on an image to reconstruct
the surface normals from the image.
# Details
Under the assumption that the object is illuminated from directly above ([0,0,1])
and that the surface is Lambertian, the eikonal equation, shown below, can be used
to describe the surface in an image.
```math
E_{i,j}=\frac{E_{max}}{\sqrt{p_{i,j}^2+q_{i,j}^2+1}}
```
As the surface must be integrable, the integrability condition; ``p_y=q_x`` must
hold on the surface. Using this a discreate itegrability function can be derived
in the form:
```math
p_{i,j+1}-p_{i,j}=q_{i+1,j}-q_{i,j}
```
A discreate functional for the energy can now be derived from these in the form:
```math
\begin{aligned}
\epsilon&=\sum_{i,j\in\Omega}\left[\frac{E_{max}}{\sqrt{p_{i,j}^2+q_{i,j}^2+1}}-E_{i,j}\right]^2\\
&+\lambda_i\sum_{i,j\in\Omega}[(p_{i,j+1}-p_{i,j})-(q_{i+1,j}-q_{i,j})]^2
\end{aligned}
```
Finally a smoothness term can be introduced as below. This both aids in the quality
of the reconstruction and reduces the number of minima in the functional helping
to ensure a good solution is found.
```math
\begin{aligned}
\epsilon&=\sum_{i,j\in\Omega}\left[\frac{E_{max}}{\sqrt{p_{i,j}^2+q_{i,j}^2+1}}-E_{i,j}\right]^2\\
&+\lambda_i\sum_{i,j\in\Omega}[(p_{i,j+1}-p_{i,j})-(q_{i+1,j}-q_{i,j})]^2\\
&+\lambda_i\sum_{i,j\in\Omega}[(p_{i,j+1}-p_{i,j})^2+(p_{i+1,j}-p_{i,j})^2+(q_{i+1,j}-q_{i,j})^2+(q_{i,j+1}-q_{i,j})^2]
\end{aligned}
```
The problem is then minimized using gradient decsent method to attampt to find
a solution minimizing the energy function. Gradient descent is terminated when
``|\nabla\epsilon|<\beta\sqrt{2|\Omega|}``.
# Parameters
## `p`
An `AbstractArray` which contains initial solution for p, the x gradient of the
surface. Better inital solution will produce better results.
## `q`
An `AbstractArray` which contains initial solution for q, the y gradient of the
surface. Better inital solution will produce better results.
## `integration_factor`
A `Real` which controls the impact of the integration term on the energy function.
Can be used to tune the final solution. Defults to 10.
## `smoothness`
A `Real` which controls the impact of the smoothness term on the energy function.
Can be used to tune the final solution. Defults to 20.
## `β`
A `Real` which controls the value of ``|\nabla\epsilon|`` which will terminate
the algorithm. If algorithm terminates on first iteration value will need reducing
and if algorithm fails to converge the value will need increasing. Defults to 1.0.
## `integration_scheme`
A `AbstractIntegrator` which will be used to integrate the calculated gradients.
Defults to `Horn()`.
# Example
The following example demonstrates the use of the `Deterministic` method.
```julia
using ShapeFromShading, Makie, Images

#generate synthetic image and initial solution
img = generate_surface(SynthSphere(38))
pin, qin = synthetic_gradient(SynthSphere(38))

#calculate the heightmap
deterministic = Deterministic(p=pin,q=qin, integration_factor = 10, smoothness = 50, β = 1.0, integration_scheme = Horn())
Z,p,q = deterministic(img)

#normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(abs.(Z))

#display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
surface(r, r, Z)
```
# Reference
[1] A. Crouzil, X. Descombes and J. Durou, "A multiresolution approach for shape from shading coupling deterministic and stochastic optimization," in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 25, no. 11, pp. 1416-1421, Nov. 2003. [doi: 10.1109/TPAMI.2003.1240116](https://doi.org/10.1109/TPAMI.2003.1240116 )
"""
function (algorithm::Deterministic)(img::AbstractArray)
    # Setup initial variables
    E = Float64.(img)
    x₀ = zeros(Float64,first(size(img)),last(size(img)),2)
    x₀[:,:,1] .= algorithm.p
    x₀[:,:,2] .= algorithm.q
    λᵢ = algorithm.integration_factor
    λₛ = algorithm.smoothness
    β = algorithm.β

    # Calcualte derivatives
    od = OnceDifferentiable(pq -> f(E, pq, λᵢ, λₛ), x₀; autodiff = :forward);
    a = zeros(size(x₀))

    # Run gradient descent
    res = optimize(od, x₀, GradientDescent(), Optim.Options(show_trace = true, g_tol = β*sqrt(2*length(E))))

    # Get minimised gradients
    out = Optim.minimizer(res)
    p = out[:,:,1]
    q = out[:,:,2]
    Z = algorithm.integration_scheme(p,q)
    return Z,p,q
end
