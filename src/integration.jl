@doc raw"""
```
Frankot()
```
Defines the Frankot integrator which contians the Frankot-Chellappa method of
integration. The Frankot-Chellappa method is a fast, reliable method for
integration surface normals while enforcing integrability using Fourier methods.
# Output
`Frankot()` returns a Frankot integrator which can then be called to run the
Frankot-Chellappa method on a gradient field.
# Details
Frankot-Chellappa method uses Fourier methods to attempt to solve the Poission
equation ``\nabla^2z = \partial_up + \partial_vq``. By taking the Fourier transform
of both sides we get:
```math
−(\omega^2_u + \omega^2_v)\hat{z}(\omega_u, \omega_v) = \imath \omega_u\hat{p}
(\omega_u, \omega_v) + \imath \omega_v\hat{q}(\omega_u, \omega_v)
```
By rearranging the above equation we arrive at an equation for ``\hat{z}``;
```math
\hat{z}(\omega_u, \omega_v) = \frac{\omega_u\hat{p}(\omega_u, \omega_v) +
\omega_v\hat{q}(\omega_u, \omega_v)}{\imath(\omega^2_u + \omega^2_v)}
```
From which the final surface can be found by taking the inverse Fourier transform
of ``\hat{z}``.

Due to the way ``(\omega_u, \omega_v)`` is defined the algorithm works best when
the input dimensions are odd in length. To accommodate this the integrator will
pad the edge of the inputs if they are even before running the algorithm. This
padding will be removed before returning a value hence output size will be
unaffected.
# Parameters
`Frankot` integrator take no parameters.
# Example
The following example demonstrates  the use of the `Frankot` integrator.
```julia
using ShapeFromShading, Makie

# Generate synthetic gradients
p, q = synthetic_gradient(SynthSphere(38), img_size = 151)

# Create a Frankot() integrator
frankot = Frankot()

# Calculate the heightmap from the gradients
Z = frankot(p, q)

# Normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(Z)

# Display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
surface(r, r, Z)
```
# Reference
[1] R. T. Frankot and R. Chellappa, "A method for enforcing integrability in shape from shading algorithms," in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 10, no. 4, pp. 439-451, July 1988. [doi: 10.1109/34.3909](https://doi.org/10.1109/34.3909 )
"""
function (scheme::Frankot)(pin::AbstractArray, qin::AbstractArray)
    # Resize grid to be odd dimensions. (Requared to setup transofrm values)
    M, N = size(pin)
    p = zeros(Complex{Float64}, M + 1 - (M%2), N + 1 - (N%2))
    p[1:M, 1:N] = Complex{Float64}.(pin)
    q = zeros(Complex{Float64}, M + 1 - (M%2), N + 1 - (N%2))
    q[1:M, 1:N] = Complex{Float64}.(qin)

    wx, wy = setup_transform_values(N,M)
    fft!(p)
    fft!(q)
    Z = zeros(Complex{Float64}, size(p))
    Z = (-1im .* wx .* p .+ 1im .* wy .* q) ./ (2 .* π .* (wx.^2 .+ wy.^2 .+ eps()))
    ifft!(Z)
    Z = abs.(Z)
    return Z[1:M,1:N]
end

@doc raw"""
```
Path()
```
Creates a `Path()` integrator which utilises the average of two path integrals
along varying paths. Each path integral reconstructs the surface with
accumulating error along the path, hence averaging two different paths can
minimise this error, although the method still suffers if the gradient field is
not integrable at some points.
# Output
`Path()` returns a Path integrator which can then be called to integrate a
gradient field.
# Details
Under the assumption that the surface normals are approximately integrable
everywhere (``\frac{\partial p}{\partial y}\approx\frac{\partial q}{\partial x}``),
then surface can be reconstructed using the path integral defined as:
```math
z(x,y)=\oint_c\left(\frac{\partial z}{\partial x},\frac{\partial z}{\partial y}\right)\cdot dl
```
Which can be broken into two integrals representing the value at each point on
the surface as shown below for a path which integrates along the first column
then along the row.
```math
z(u,v)=\int_0^v\frac{\partial z}{\partial y}(0,y)dy + \int_0^u\frac{\partial z}{\partial x}(x,v)dx
```
The second path used in the algorithm is simply the transpose of the first,
integrating along the first row then down the column represented mathematically as:
```math
z(u,v)=\int_0^u\frac{\partial z}{\partial x}(x,0)dx + \int_0^v\frac{\partial z}{\partial y}(u,y)dy
```
The algorithm can be written, then discretised as shown below:
```math
\begin{gathered}
z(u,v)=\frac{1}{2}\left(\int_0^v\frac{\partial z}{\partial y}(0,y)dy + \int_0^u\frac{\partial z}{\partial x}(x,v)dx + \int_0^u\frac{\partial z}{\partial x}(x,0)dx + \int_0^v\frac{\partial z}{\partial y}(u,y)dy\right)\\
z(u,v)=\frac{1}{2}\left(\sum_{i=0}^vq(0,i) + \sum_{j=0}^up(j,v) + \sum_{j=0}^up(j,0) + \sum_{i=0}^vq(u,i)\right)\\
z(u,v)=\frac{1}{2}\left(\sum_{i=0}^v(q(0,i) + q(u,i)) + \sum_{j=0}^u(p(j,0) + p(j,v))\right)\\
\end{gathered}
```
It is important to note as mentioned above if there are non-integrable points in
the normal field then artefacts can appear in the reconstruction. This is
seen in the example below where the otherwise smooth sphere appears "spiky".
This can be corrected post reconstruction by smoothing but ideally a different
integrator should be used.
# Parameters
`Path` integrator take no parameters.
# Example
The following example demonstrates the use of the `Path` integrator.
```julia
using ShapeFromShading, Makie

# Generate synthetic gradients
p, q = synthetic_gradient(SynthSphere(38), img_size = 151)

# Create a Path() integrator
path = Path()

# Calculate the heightmap from the gradients
Z = path(p, q)

# Normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(Z)

# Display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
surface(r, r, Z)
```
# Reference
[1] D. Forsyth and J. Ponce, Computer vision: a modern approach. Upper Saddle River, N.J: Prentice Hall, 2003, pp. 84-86.
"""
function (scheme::Path)(p::AbstractArray, q::AbstractArray)
    R, C = size(p)
    Z₁ = zeros(Float64, size(p))
    Z₂ = zeros(Float64, size(p))

    #path 1
    for i = 2:R
        Z₁[i,1] = Z₁[i-1,1] - q[i,1]
    end
    for i = 1:R
        for j = 2:C
            Z₁[i,j] = Z₁[i,j-1] + p[i,j]
        end
    end

    #path 2
    for i = 2:C
        Z₂[1,i] = Z₂[1,i-1] + p[1,i]
    end
    for i = 2:R
        for j = 1:C
            Z₂[i,j] = Z₂[i-1,j] - q[i,j]
        end
    end
    Z = (Z₁ .+ Z₂) / 2
    return Z
end

@doc raw"""
```
SplitPath()
```
Creates a `SplitPath()` integrator which utilizes the average of two path integrals
along varying paths averaging the value at each step. Each path integral
reconstructs the surface with accumlating error along the path, hence averaging
two different paths at each step reduces the global error at the cost of local
error, although the method still suffers if the gradient field is not integrable
at some points it does less so the `Path()` from which it extends.
# Output
`SplitPath()` returns a SplitPath integrator which can then be called to integrate
a gradient field.
# Details
Under the assumption that the surface normals are approximately integrable
everywhere (``\frac{\partial p}{\partial y}\approx\frac{\partial q}{\partial x}``),
then surface can be reconstructed using the path integral defined as:
```math
z(x,y)=\oint_c\left(\frac{\partial z}{\partial x},\frac{\partial z}{\partial y}\right)\cdot dl
```
By expanding on this principle and the discreate summation from `Path()` we can
arrive at the discreate expresion for the value at each point, assuming all
values prior to that point have been calculated, as follows:
```math
z_{u,v} = \frac{1}{2}(z_{u-1,v}+p_{u-1,v}+z_{u,v-1}+q_{u,v-1})
```
As with other similar methods (see `Horn()`) care must be taken with regards to
boundaries which can be calculated, to a constant value ``z(0,0)`` which is assumed
to be the zero point, using:
```math
\begin{gathered}
z_{u,0} = z_{u-1,0}+p_{u-1,0}\\
z_{0,v} = z_{0,v-1}+q_{0,v-1}
\end{gathered}
```
It is important to note as mentioned above if there are non-integrable points in
the normal field then artefacts can appear in the reconstruction. These errors
gradually average out but will lead to "streaks" appearing in the reconstruction.
This is seen in the example below where the otherwise smooth sphere appears has
ripple like structures pointing toward to top right corner. This can be corrected
post reconstruction by smoothing but ideally a different integrator should be used.
It is also interesting to note the parallels between this method and the Horn
and Brooks method, with this method being effectively the forward component of
Horn's method. As such this algorithm provided a useful middle ground between
direct integration algorithms and iterative algorithms such as the Horn and Brooks
method.
# Parameters
`SplitPath` integrator take no parameters.
# Example
The following example demonstrates the use of the `SplitPath` integrator.
```julia
using ShapeFromShading, Makie

# Generate synthetic gradients
p, q = synthetic_gradient(SynthSphere(38), img_size = 151)

# Create a Path() integrator
splitPath = SplitPath()

# Calculate the heightmap from the gradients
Z = splitPath(p, q)

# Normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(Z)

# Display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
surface(r, r, Z)
```
# Reference
[1] D. Forsyth and J. Ponce, Computer vision: a modern approach. Upper Saddle River, N.J: Prentice Hall, 2003, pp. 84-86.
[2] B. Horn and M. Brooks, "The variational approach to shape from shading", Computer Vision, Graphics, and Image Processing, vol. 33, no. 2, pp. 174-208, 1986. [doi: 10.1016/0734-189x(86)90114-3](https://doi.org/10.1016/0734-189x(86)90114-3 )
"""
function (scheme::SplitPath)(p::AbstractArray, q::AbstractArray)
    R, C = size(p)
    Z = zeros(Float64, size(p))
    for i = 2:R
        Z[i,1]=Z[i-1,1] + q[i,1]
        Z[1,i]=Z[1,i-1] - p[1,i]
    end
    for i = 2:R
        for j = i:R
            Z[i,j] = ((Z[i,j-1] + p[i,j]) + (Z[i-1,j] - q[i,j])) / 2
            Z[j,i] = ((Z[j,i-1] + p[j,i]) + (Z[j-1,i] - q[j,i])) / 2
        end
    end
    return Z
end

@doc raw"""
```
Horn(ϵ::Real = 1.0, max_iter::Real = 10000)
```
Implements the Horn and Brook's method of integrating surface normals. This
algorithm offers an iterative solution to the Poisson equation describing the
surface providing good reconstructions under most conditions.
# Output
`Horn()` returns a Horn integrator which can then be called to integrate
a gradient field.
# Details
The Horn and Brook's method attempts to solve the Poisson equation ``\nabla^2z = \partial_up + \partial_vq``
by the discretization given below:
```math
z_{u+1,v}+z_{u,v+1}+z_{u-1,v}+z_{u,v-1}-4z_{u,v}=\frac{p_{u+1,v}-p_{u-1,v}}{2}+\frac{q_{u,v+1}-q_{u,v-1}}{2}
```
Which can be rearranged to give the iterative scheme provided by:
```math
z_{u,v}^{k+1}= \frac{z_{u+1,v}^k + z_{u,v+1}^k + z_{u-1,v}^k + z_{u,v-1}^k}{4} - \frac{p_{u+1,v}-p_{u-1,v}}{8} - \frac{q_{u,v+1}-q_{u,v-1}}{8}
```
This scheme will always converge to a solution however the rate of convergence
may depend upon the initial solution. This implementation will initilize with a
zero solution. Neumann boundary conditions are imposed at the edges where the
scheme would otherwise go out of bounds.
# Parameters
The function parameters are described in more detail below.
##  `Max_iter`:
An `Int` which controls the number of iterations the algorithm will run for.
## `ϵ`:
A `Real` representing the distance between pixels. This will Control how tall
the final reconstruction is relative the array grid.
# Example
The following example demonstrates the use of the `Horn` integrator.
```julia
using ShapeFromShading, Makie

# Generate synthetic gradients
p, q = synthetic_gradient(SynthSphere(38), img_size = 151)

# Create a Horn() integrator
horn = Horn(ϵ = 0.03, max_iter = 10000)

# Calculate the heightmap from the gradients
Z = horn(p, q)

# Display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
surface(r, r, Z)
```
# Reference
[1] B. Horn and M. Brooks, "The variational approach to shape from shading", Computer Vision, Graphics, and Image Processing, vol. 33, no. 2, pp. 174-208, 1986. [doi: 10.1016/0734-189x(86)90114-3](https://doi.org/10.1016/0734-189x(86)90114-3 )
"""
function (scheme::Horn)(p::AbstractArray, q::AbstractArray)
    iter = copy(scheme.max_iter)
    ϵ = copy(scheme.ϵ)
    Z = zeros(Float64, size(p))
    Zᵏ⁺¹ = zeros(Float64, size(p))
    h = zeros(Float64, size(p))
    v = zeros(Float64, size(p))
    hv = zeros(Float64, size(p))
    R, C = size(p)
    for i = 1:R
        Z[i,1] = p[i,1]
        Z[i,R] = p[i,R]
    end
    for i = 1:C
        Z[1,i] = q[1,i]
        Z[C,i] = q[C,i]
    end
    for i = 2:(R-1)
        for j = 2:(C-1)
            h[i,j] = (p[i,j+1] - p[i,j-1]) / 2.0
            v[i,j] = (q[i+1,j] - q[i-1,j]) / 2.0
            hv[i,j] = (ϵ / 4) * (h[i,j] - v[i,j])
        end
    end
    for k = 1:iter
        copyto!(Zᵏ⁺¹, Z)
        for i = 2:(R-1)
            for j = 2:(C-1)
                Zᵏ⁺¹[i,j] = (Z[i-1,j] + Z[i+1,j] + Z[i,j-1] + Z[i,j+1]) / 4.0
                Zᵏ⁺¹[i,j] = Zᵏ⁺¹[i,j] - hv[i,j]
            end
        end
        Z = Zᵏ⁺¹
    end
    return Z
end

@doc raw"""
```
Durou(ϵ::Real = 1.0, max_iter::Real = 1000)
```
Implements the Durou and Courteille method of integrating surface normals. This
algorithm offers an iterative solution to the Poisson equation describing the
surface extending Horn and Brook's method by improving the boundary approximation
and providing good reconstructions under most conditions.
# Output
`Durou()` returns a Durou integrator which can then be called to integrate
a gradient field.
# Details
TheDurou and Courteille's method attempts to solve the Poisson equation
``\nabla^2z = \partial_up + \partial_vq`` by the discretization given below:
```math
z_{u+1,v}+z_{u,v+1}-2z_{u,v}=\frac{p_{u+1,v}+p_{u,v}}{2}+\frac{q_{u,v+1}+q_{u,v}}{2}
```
Which can be rearranged to give the iterative scheme provided by:
```math
z_{u,v}^{k+1}= \frac{z_{u+1,v}^k + z_{u,v+1}^k}{2} - \frac{p_{u+1,v}+p_{u,v}}{4} - \frac{q_{u,v+1}+q_{u,v}}{4}
```
This scheme will always converge to a solution however the rate of convergence
may depend upon the initial solution. This implementation will initialize with a
zero solution. Natural boundary conditions are imposed at the edges using the
condition ``\partial_uz-p+\partial_v-q=0``. Although faster then the Horn and Brook's
method and better at handling boundaries, it can generate a worse solution under
some conditions.
# Parameters
The function parameters are described in more detail below.
##  `Max_iter`:
An `Int` which controls the number of iterations the algorithm will run for.
the range [0,1].
## `ϵ`:
A `Real` representing the distance between pixels. This will Control how tall
the final reconstruction is relative the array grid.
# Example
The following example demonstrates the use of the `Durou` integrator.
```julia
using ShapeFromShading, Makie

# Generate synthetic gradients
p, q = synthetic_gradient(SynthSphere(38), img_size = 151)

# Create a Durou() integrator
durou = Durou(ϵ = 0.03, max_iter = 10000)

# Calculate the heightmap from the gradients
Z = durou(p, q)

# Display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
surface(r, r, Z)
```
# Reference
[1] Y. Quéau, J. Durou and J. Aujol, "Normal Integration: A Survey", Journal of Mathematical Imaging and Vision, vol. 60, no. 4, pp. 576-593, 2017. [doi: 10.1007/s10851-017-0773-x](https://doi.org/10.1007/s10851-017-0773-x )
"""
function (scheme::Durou)(p::AbstractArray, q::AbstractArray)
    iter = copy(scheme.max_iter)
    ϵ = copy(scheme.ϵ)
    Z = zeros(Float64, size(p))
    Zᵏ⁺¹ = zeros(Float64, size(p))
    h = zeros(Float64, size(p))
    v = zeros(Float64, size(p))
    hv = zeros(Float64, size(p))
    R, C = size(p)
    for i = 1:(R-1)
        for j = 1:(C-1)
            h[i,j] = (p[i,j+1] + p[i,j]) / 2
            v[i,j] = (q[i+1,j] + q[i,j]) / 2
            hv[i,j] = (ϵ / 2) * (h[i,j] - v[i,j])
        end
    end
    for k = 1:iter
        copyto!(Zᵏ⁺¹, Z)
        for i = 1:(R-1)
            for j = 1:(C-1)
                Zᵏ⁺¹[i,j] = (Z[i+1,j] + Z[i,j+1]) / 2
                Zᵏ⁺¹[i,j] = Zᵏ⁺¹[i,j] - hv[i,j]
            end
        end
        Z = Zᵏ⁺¹
    end
    return Z
end

# Generates matricies Dᵤ⁺, Dᵤ⁻, Dᵥ⁺ and Dᵥ⁻ as per eq 11, 12, 13, 14 and 17
function gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻,T,mask)
    U,V = size(mask)
    for i in CartesianIndices(mask)
        #u+
        if i[1] + 1 <= U && mask[i] == 1 && mask[i[1]+1,i[2]] == 1
            Dᵤ⁺[i[1]+(i[2]-1)*U, i[1]+(i[2]-1)*U] = -1.0
            Dᵤ⁺[i[1]+(i[2]-1)*U, i[1]+(i[2]-1)*U+1] = 1.0
        end
        #u-
        if i[1] - 1 > 0 && mask[i] == 1 && mask[i[1]-1,i[2]] == 1
            Dᵤ⁻[i[1]+(i[2]-1)*U, i[1]+(i[2]-1)*U] = 1.0
            Dᵤ⁻[i[1]+(i[2]-1)*U, i[1]+(i[2]-1)*U-1] = -1.0
        end
        #v+
        if i[2] + 1 <= V && mask[i] == 1 && mask[i[1],i[2]+1] == 1
            Dᵥ⁺[i[1]+(i[2]-1)*V, i[1]+(i[2]-1)*V] = -1.0
            Dᵥ⁺[i[1]+(i[2]-1)*V, i[1]+(i[2])*V] = 1.0
        end
        #v-
        if i[2] - 1 > 0 && mask[i] == 1 && mask[i[1],i[2]-1] == 1
            Dᵥ⁻[i[1]+(i[2]-1)*V, i[1]+(i[2]-1)*V] = 1.0
            Dᵥ⁻[i[1]+(i[2]-1)*V, i[1]+(i[2]-2)*V] = -1.0
        end
    end
end

@doc raw"""
```
Quadratic(z::AbstractArray, λ::AbstractArray = fill(10.0^-6, size(z)), mask::AbstractArray = fill(1.0, size(z)))
```
Implements the quadratic variational least squared method proposed by Aujol, Durou
and Quéau. The algorithm solves the least squares system generated from minimizing
a fidelity term ``\mathcal{F}(z) = \iint_{(u,v)\in \Omega}\Phi(||\nabla z(u,v) - g(u,v)||)dudv``
and a regularization term ``\mathcal{R}(z) = \iint_{(u,v)\in \Omega}\lambda \left[z(u,v)-z^0(u,v)\right]^2dudv``
where ``\Phi(s) = s^2``. This method is able to quickly produce good solutions on
a smooth surface and can easily handle non-rectangular domains or sub-divided
domains and can provide a good starting solution for other algorithms.
# Output
`Quadratic()` returns a Quadratic integrator which can then be called to run the
quadratic method on a gradient field.
# Details
As mentioned above the quadratic variational least squared method aims to minimize
a fidelity term ``\mathcal{F}(z) = \iint_{(u,v)\in \Omega}\Phi(||\nabla z(u,v) - g(u,v)||)dudv``
and a regulization term ``\mathcal{R}(z) = \iint_{(u,v)\in \Omega}\lambda \left[z(u,v)-z^0(u,v)\right]^2dudv``
where ``\Phi(s) = s^2``.

This leads to the minimization problem given by:
```math
min\iint_{(u,v)\in \Omega}||\nabla z(u,v)-\mathbf{g}(u,v)||^2+\lambda(u,v)\left[z(u,v)-z^0(u,v)\right]^2dudv
```
By discretizing the problem we can arrive at the functional:
```math
\begin{aligned}
E(z) = \frac{1}{2}\Big(\sum_{(u,v)\in \Omega_u^+}[\partial_u^+z_{u,v}-p_{u,v}]^2 \\+ \sum_{(u,v)\in \Omega_u^-}[\partial_u^-z_{u,v}-p_{u,v}]^2 \\+ \sum_{(u,v)\in \Omega_v^+}[\partial_v^+z_{u,v}-q_{u,v}]^2 \\+ \sum_{(u,v)\in \Omega_v^-}[\partial_v^-z_{u,v}-q_{u,v}]^2\Big) \\+ \sum_{(u,v)\in \Omega}[z_{u,v}-z^0_{u,v}]^2
\end{aligned}
```
where ``\Omega_u^+`` represents the domain where ``(u,v)\in\{(u,v)\in\Omega|(u,v)(u+1,v)\in\Omega\}`` etc..
Using this definition the discrete differences can be converted into matrix form,
where ``p_{u,v}=z_{u+1,v}-z_{u,v}`` is the forward difference in the u direction etc.
This data is then stacked into three vectors; ``\mathbf{z},\mathbf{p},\mathbf{q}\in\R^{|\Omega|}``.
Thus the matrix reresenting each of these is defined as below where m(i) is the
mapping of the ith element of this vector to its corresponding point in ``(u,v)``
and ``D_u^+`` is a ``|\Omega|\times|\Omega|`` matrix.
```math
D_u^+[i,j]=\begin{cases}
   0 &\text{if } m(i)\notin\Omega_u^+ \text{or } j \ne i \text{or } j \ne i+1\\
   -1 &\text{if } j = i \text {and } m(i)\in\Omega_u^+\\
   1 &\text{if } j = i+1 \text {and } m(i)\in\Omega_u^+
\end{cases}
```
For a 2X2 domain this looks like:
```math
D_u^+ = \begin{bmatrix}
   -1 & 1 & 0 & 0 \\
   0 & 0 & 0 & 0 \\
   0 & 0 & -1 & 1 \\
   0 & 0 & 0 & 0
\end{bmatrix}
```
The other three discrete differences matrices are similarly defined from there
definitions to be ``D_u^-``, ``D_v^+`` and ``D_v^-``. These can then be used to
redefine the minimization problem to be in the form:
```math
E(\mathbf{z})=\frac{1}{2}\left(||D_u^+\mathbf{z}-\mathbf{p}||^2+||D_u^-\mathbf{z}-\mathbf{p}||^2+||D_v^+\mathbf{z}-\mathbf{q}||^2+||D_v^-\mathbf{z}-\mathbf{q}||^2\right)+||\Lambda(\mathbf{z}-\mathbf{z}^0)||^2
```
Where ``\Lambda`` is the ``|\Omega|\times|\Omega|`` diagonal matrix containing
the values of ``\sqrt{\lambda_{u,v}}``. Using the above definitions the negative
Laplacian matrix can then be defined as:
```math
L=\frac{1}{2}[D_u^{+\top}D_u^++D_u^{-\top}D_u^-+D_v^{+\top}D_v^++D_v^{-\top}D_v^-]
```
Finally the least minimization problem can be represented in the form of a least
squares problem of the form ``A\mathbf{z}=\mathbf{b}`` where:
```math
\begin{gathered}
A=L+\Lambda^2\\
\mathbf{b}=\frac{1}{2}\left[D_u^{+\top}+D_u^{-\top}\right]\mathbf{p}+\frac{1}{2}\left[D_v^{+\top}+D_v^{-\top}\right]\mathbf{q}+\Lambda^2\mathbf{z}^0\\
=D_u\mathbf{p}+D_v\mathbf{q}+\Lambda^2\mathbf{z}^0
\end{gathered}
```
This system is then solved using a standard conjugate gradient algorithm where
the initialization has only a slight impact on the runtime and no impact on the
final solution. The algorithm provides good results on smooth surfaces but
struggles in the presence of discontinuities.
# Parameters
## `z`:
An `AbstractArray` which defines the value of ``z^0`` the initial solution and
prior to be used in the regularization term. Must be provided.
## `λ`:
An `AbstractArray` the same size as `z`, defaulting to ``10.0^{-6}`` everywhere.
This defines theregulization weight at each point. Large values will force the
algorithm to keep the solution near to ``z^0`` at that position. Can be used to
keep the solution near the initial solution or guide the solution to a certain
known value at points (i.e. known maxima and minima). This value should be set
uniformly small otherwise.
## `mask`:
An `AbstractArray` the same size as `z`, which guides the algorithm as to where
the valid domain is. Values of `1` will be in the domain ``\Omega`` while other
values will be ignored and set to ``z^0``. This can be used to
integrate over sub-domain or to segment the domain into parts. The `gen_mask()`
funtion can be used to generate a mask which will remove non-integrable regions
dramatically improving the solution under most condition at the cost of not
integrating the entire solution.
# Example
The following example demonstrates the use of the `Quadratic` integrator.
```julia
using ShapeFromShading, Makie

# Generate synthetic gradients
p, q = synthetic_gradient(Prism(75), img_size = 151)

# Create a Quadratic() integrator
quadratic = Quadratic(z=zeros(size(p)))
quadraticMasked = Quadratic(z=zeros(size(p)), mask=gen_mask(p,q,1.0)[:,:,1])

# Calculate the heightmap from the gradients
Z = quadratic(p, q)
Z2 = quadraticMasked(p, q)

# Normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(Z)
Z2 = Z2./maximum(Z2)

# Display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
vbox(surface(r, r, Z), surface(r, r, Z2))
```
# Reference
[1] Y. Quéau, J. Durou and J. Aujol, "Variational Methods for Normal Integration", Journal of Mathematical Imaging and Vision, vol. 60, no. 4, pp. 609-632, 2017. [doi: 10.1007/s10851-017-0777-6](https://doi.org/10.1007/s10851-017-0777-6 )
"""
function (scheme::Quadratic)(pIn::AbstractArray, qIn::AbstractArray)
    z⁰ = copy(scheme.z)
    λ = copy(scheme.λ)
    mask = copy(scheme.mask)
    z = copy(z⁰)
    mask = rotr90(mask)
    λ = rotr90(λ)
    p = rotr90(copy(pIn)).*mask
    q = rotr90(copy(qIn)).*mask
    T = first(size(z))
    index = CartesianIndices(z)
    index = reshape(index, length(index),1)

    Dᵤ⁺ = spzeros(Float64, T^2, T^2)
    Dᵤ⁻ = spzeros(Float64, T^2, T^2)
    Dᵥ⁺ = spzeros(Float64, T^2, T^2)
    Dᵥ⁻ = spzeros(Float64, T^2, T^2)
    gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻, T, mask)

    # Calculate Dᵤ, Dᵥ and L as per eq 23 and 24
    Dᵤ = 0.5*(transpose(Dᵤ⁺) .+ transpose(Dᵤ⁻))
    Dᵥ = 0.5*(transpose(Dᵥ⁺) .+ transpose(Dᵥ⁻))
    L = 0.5*((transpose(Dᵤ⁺) * Dᵤ⁺) .+ (transpose(Dᵤ⁻) * Dᵤ⁻) .+ (transpose(Dᵥ⁺) * Dᵥ⁺) .+ (transpose(Dᵥ⁻) * Dᵥ⁻))

    Λ = spzeros(Float64,T^2,T^2)
    A = spzeros(Float64,T^2,T^2)
    P = zeros(Float64,T^2)
    Q = zeros(Float64,T^2)
    Z⁰ = zeros(Float64,T^2)
    Z = zeros(Float64,T^2)

    # Vectorize inputs
    for j in CartesianIndices(z)
        Λ[j[1]+(j[2]-1)*T,j[1]+(j[2]-1)*T] = sqrt(λ[j])
        P[j[1]+(j[2]-1)*T] = p[j]
        Q[j[1]+(j[2]-1)*T] = q[j]
        Z⁰[j[1]+(j[2]-1)*T] = z⁰[j]
        Z[j[1]+(j[2]-1)*T] = z[j]
    end

    # Calculate A and b as per eq 23 abd 24
    A = L .+ (Λ^2)
    b = Dᵤ*P .+ Dᵥ*Q .+ (Λ^2)*Z⁰

    # Solve system
    pl = AMGPreconditioner{RugeStuben}(A)
    cg!(Z, A, b, tol=10.0^-4, Pl = pl)

    for j in eachindex(index)
        z[index[j]] = Z[j]
    end
    z = z.*mask
    z = rotl90(z)
    return z
end

@doc raw"""
```
TotalVariation(z::AbstractArray, α::Real = 1.0, λ::AbstractArray = fill(10.0^-6, size(z)), mask::AbstractArray = fill(1.0, size(z)), max_iter::Int = 100)
```
Implements the total variational method proposed by Aujol, Durou and Quéau. The
algorithm solves the same minimization problem as th `Quadratic`method exept the
fidelity terms function ``\Phi(s)=||s||_{L_1}``. This method is able to produce
good solutions on a smooth and piecewise smooth surface and can easily handle
non-rectangular domains or sub-divided domains.
# Output
`TotalVariation()` returns a TotalVariation integrator which can then be called to run the
quadratic method on a gradient field.
# Details
As discussed above this algorithm used the same fidelity and regularization terms
as the `Quadratic` method exept the ``\Phi(s)=s^2`` term is replaced with
``\Phi(s)=||s||_{L_1}``. This leads to the minimization problem defined by:
```math
min\iint_{(u,v)\in \Omega}||\nabla z(u,v)-\mathbf{g}(u,v)||+\lambda(u,v)\left[z(u,v)-z^0(u,v)\right]^2dudv
```
By considering the four posible discreatisations of ``∇z(u,v)`` we can generate
four posible domains to consider given by; ``\Omega^{UV}=\Omega_u^U\cup\Omega_v^V, (U,V)\in\{+,-\}^2``
where ``\{+,-\}^2`` refers to all posible combinations of +,-.
Using these the following discreate functional can be generated:
```math
\begin{aligned}
E(\mathbf{z})=\frac{1}{4}\Big(&\sum_{(u,v)\in\Omega^{++}}\sqrt{[\partial_u^+z_{u,v}-p_{u,v}]^2+[\partial_v^+z_{u,v}-q_{u,v}]^2}\\
+&\sum_{(u,v)\in\Omega^{+-}}\sqrt{[\partial_u^+z_{u,v}-p_{u,v}]^2+[\partial_v^-z_{u,v}-q_{u,v}]^2}\\
+&\sum_{(u,v)\in\Omega^{-+}}\sqrt{[\partial_u^-z_{u,v}-p_{u,v}]^2+[\partial_v^+z_{u,v}-q_{u,v}]^2}\\
+&\sum_{(u,v)\in\Omega^{--}}\sqrt{[\partial_u^-z_{u,v}-p_{u,v}]^2+[\partial_v^-z_{u,v}-q_{u,v}]^2}\Big)\\
+&\sum_{(u,v)\in\Omega}\lambda_{u,v}\left[z_{u,v}-z^0_{u,v}\right]^2
\end{aligned}
```
Which simplifies to the minimization problem:
```math
\begin{gathered}
min\frac{1}{4}\sum_{(U,V)\in\{+,-\}^2}\sum_{(u,v)\in\Omega^{UV}}||\mathbf{r}_{(u,v)}^{UV}||+\sum_{(u,v)\in\Omega}\lambda_{u,v}\left[z_{u,v}-z^0_{u,v}\right]^2\\
\mathbf{r}_(u,v)^{UV}=\nabla^{UV}z_{u,v}-\mathbf{g}_{u,v}
\end{gathered}
```
This leads to the optimization scheme using an ADMM algorithm defined by:
```math
\begin{aligned}
z^{(k+1)}=&min\frac{\alpha}{8}\sum_{(U,V)\in\{+,-\}^2}\sum_{(u,v)\in\Omega^{UV}}||\nabla^{UV}z_{u,v}-(\mathbf{g}_{u,v}+\mathbf{r}_{(u,v)}^{UV^{(k)}}-\mathbf{b}_{(u,v)}^{UV^{(k)}})||\\&+\sum_{(u,v)\in\Omega}\lambda_{u,v}\left[z_{u,v}-z^0_{u,v}\right]^2\\
\mathbf{r}_{(u,v)}^{UV^{(k+1)}}=&min\frac{\alpha}{8}||\mathbf{r}-(\nabla^{UV}z_{u,v}-\mathbf{g}_{u,v}+\mathbf{b}_{(u,v)}^{UV^{(k)}})||+||\mathbf{r}||\\
\mathbf{b}_{(u,v)}^{UV^{(k+1)}}=&\mathbf{b}_{(u,v)}^{UV^{(k)}})+\nabla^{UV}z_{u,v}-\mathbf{g}_{u,v}-\mathbf{r}_{(u,v)}^{UV^{(k+1)}}
\end{aligned}
```
The z update then can be solved using the linear system defined below, where D_{u,v}^{U,V}
and ``\Lambda`` are the same at those defined in [`Quadratic`](@ref).
```math
\begin{aligned}
    A_{TV}&\mathbf{z}^{(k+1)}=b_{TV}^{(k)}\\
    A_{TV}&=\frac{\alpha}{8}\sum_{(U,V)\in\{+,-\}^2}\left[D_u^{U\top}D_u^{U}+D_v^{V\top}D_v^{V}\right] + \Lambda^2\\
    b_{TV}^{(k)}&=\frac{\alpha}{8}\sum_{(U,V)\in\{+,-\}^2}\left[D_u^{U\top}\mathbf{P}^{UV^{(k)}} + D_v^{V\top}\mathbf{Q}^{UV^{(k)}}\right]+ \Lambda^2\mathbf{z}^0\\
\end{aligned}
```
Where ``\mathbf{P}^{UV^{(k)}}, \mathbf{Q}^{UV^{(k)}}`` are the u and v components
of ``\mathbf{g}+\mathbf{r}^{UV^{(k)}}-\mathbf{b}^{UV^{(k)}}``. This can be solved
using conjugate gradient. Finally the update to ``\mathbf{r}^{UV}`` can be computed
as:
```math
\begin{aligned}
    \mathbf{r}^{UV^{(k+1)}}&=max\Big\{||\mathbf{s}_{u,v}^{UV^{(k+1)}}||-\frac{4}{\alpha},0\Big\}\frac{\mathbf{s}_{u,v}^{UV^{(k+1)}}}{||\mathbf{s}_{u,v}^{UV^{(k+1)}}||}\\
    \text{Where:}&\\
    \mathbf{s}^{UV^{(k+1)}}&=\nabla^{UV}z_{u,v}^{(k+1)}-\mathbf{g}_{u,v}+\mathbf{b}_{u,v}^{UV^{(k)}}
\end{aligned}
```
# Parameters
## `z`:
An `AbstractArray` which defines the value of ``z^0`` the initial solution and
prior to be used in the regulization term. Must be provided.
## `α`:
A `Real` with defult value of `1.0` which controls the step size. In theory this
value should have no impact on final solution but in practice larger values can
lead to worse solutions while values which are two small may lead non-convergance
in the least square update step, causing the algorithm to hang for long periods
of time. Values below `0.25` are not recomended but may work depending on domain
size and inputs.
## `λ`:
An `AbstractArray` the same size as `z` defulting to ``10.0^{-6}`` everywhere.
This defines theregulization weight at each point. Large valueas will force the
algorithm to keep the solution near to ``z^0`` at that position. Can be used to
keep the solution near the initial solution or guide the solution to a certian
known value at points (i.e. known maxima and minima). This value should be set
uniformly small otherwise.
## `mask`:
An `AbstractArray` the same size as `z`, which guides the algorithm as to where
the valid domain is. Values of `1` will be in the domain ``\Omega`` while other
values will be ignored and set to ``z^0``. This can be used to
integrate over sub-domain or to segment the domain into parts. The `gen_mask()`
funtion can be used to generate a mask which will remove non-integrable regions
dramatically improving the solution under most condition at the cost of not
integrating the entire solution.
# Example
The following example demonstrates the use of the `TotalVariation` integrator.
```julia
using ShapeFromShading, Makie

# Generate synthetic gradients
p, q = synthetic_gradient(Prism(75), img_size = 151)

# Create a TotalVariation() integrator
totalVariation = TotalVariation(z=zeros(size(p)), α=1.0)
totalVariationMasked = TotalVariation(z=zeros(size(p)), α=0.5, mask=gen_mask(p,q,1.0)[:,:,1])

# Calculate the heightmap from the gradients
Z = totalVariation(p, q)
Z2 = totalVariationMasked(p, q)

# Normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(Z)
Z2 = Z2./maximum(Z2)

# Display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
vbox(surface(r, r, Z), surface(r, r, Z2))
```
# Reference
[1] Y. Quéau, J. Durou and J. Aujol, "Variational Methods for Normal Integration", Journal of Mathematical Imaging and Vision, vol. 60, no. 4, pp. 609-632, 2017. [doi: 10.1007/s10851-017-0777-6](https://doi.org/10.1007/s10851-017-0777-6 )
"""
function (scheme::TotalVariation)(pIn::AbstractArray, qIn::AbstractArray)
    z⁰ = copy(scheme.z)
    λ = copy(scheme.λ)
    mask = copy(scheme.mask)
    α = copy(scheme.α)
    max_iter = copy(scheme.max_iter)
    z = copy(z⁰)
    mask = rotr90(mask)
    λ = rotr90(λ)
    p = rotr90(copy(pIn)).*mask
    q = rotr90(copy(qIn)).*mask
    T = first(size(z))
    index = CartesianIndices(z)
    index = reshape(index, length(index),1)

    Dᵤ⁺ = spzeros(Float64,T^2,T^2)
    Dᵤ⁻ = spzeros(Float64,T^2,T^2)
    Dᵥ⁺ = spzeros(Float64,T^2,T^2)
    Dᵥ⁻ = spzeros(Float64,T^2,T^2)
    gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻, T, mask)
    Λ = spzeros(Float64,T^2,T^2)
    A = spzeros(Float64,T^2,T^2)
    P = zeros(Float64,T^2)
    Q = zeros(Float64,T^2)
    Z⁰ = zeros(Float64,T^2)
    Z = zeros(Float64,T^2)

    # Vectorize inputs
    for j in CartesianIndices(z)
        Λ[j[1]+(j[2]-1)*T,j[1]+(j[2]-1)*T] = sqrt(λ[j])
        P[j[1]+(j[2]-1)*T] = p[j]
        Q[j[1]+(j[2]-1)*T] = q[j]
        Z⁰[j[1]+(j[2]-1)*T] = z⁰[j]
        Z[j[1]+(j[2]-1)*T] = z[j]
    end

    # Calculate A as per eq 59
    A = (α/8)*(((transpose(Dᵤ⁺) * Dᵤ⁺) .+ (transpose(Dᵥ⁺) * Dᵥ⁺)) .+ ((transpose(Dᵤ⁻) * Dᵤ⁻) .+ (transpose(Dᵥ⁻) * Dᵥ⁻)) .+ ((transpose(Dᵤ⁻) * Dᵤ⁻) .+ (transpose(Dᵥ⁺) * Dᵥ⁺)) .+ ((transpose(Dᵤ⁺) * Dᵤ⁺) .+ (transpose(Dᵥ⁻) * Dᵥ⁻))) .+ (Λ^2)
    pl = AMGPreconditioner{RugeStuben}(A)

    # Calculate b, s, r, P and Q in each combination of {+,-} as per equations 57, 61, 62
    bᵤ⁺⁺ = zeros(Float64,T^2)
    bᵥ⁺⁺ = zeros(Float64,T^2)
    sᵤ⁺⁺ = Dᵤ⁺*Z - P + bᵤ⁺⁺
    sᵥ⁺⁺ = Dᵥ⁺*Z - Q + bᵥ⁺⁺
    rᵤ⁺⁺ = max(norm(sᵤ⁺⁺) - (4/α),0)*(sᵤ⁺⁺./(norm(sᵤ⁺⁺)+eps(Float64)))
    rᵥ⁺⁺ = max(norm(sᵥ⁺⁺) - (4/α),0)*(sᵥ⁺⁺./(norm(sᵥ⁺⁺)+eps(Float64)))
    bᵤ⁺⁺ = bᵤ⁺⁺ + Dᵤ⁺*Z - P - rᵤ⁺⁺
    bᵥ⁺⁺ = bᵥ⁺⁺ + Dᵥ⁺*Z - Q - rᵥ⁺⁺
    P⁺⁺ = P + rᵤ⁺⁺ - bᵤ⁺⁺
    Q⁺⁺ = Q + rᵥ⁺⁺ - bᵥ⁺⁺

    bᵤ⁺⁻ = zeros(Float64,T^2)
    bᵥ⁺⁻ = zeros(Float64,T^2)
    sᵤ⁺⁻ = Dᵤ⁺*Z - P + bᵤ⁺⁻
    sᵥ⁺⁻ = Dᵥ⁻*Z - Q + bᵥ⁺⁻
    rᵤ⁺⁻ = max(norm(sᵤ⁺⁻) - (4/α),0)*(sᵤ⁺⁻./(norm(sᵤ⁺⁻)+eps(Float64)))
    rᵥ⁺⁻ = max(norm(sᵥ⁺⁻) - (4/α),0)*(sᵥ⁺⁻./(norm(sᵥ⁺⁻)+eps(Float64)))
    bᵤ⁺⁻ = bᵤ⁺⁻ + Dᵤ⁺*Z - P - rᵤ⁺⁻
    bᵥ⁺⁻ = bᵥ⁺⁻ + Dᵥ⁻*Z - Q - rᵥ⁺⁻
    P⁺⁻ = P + rᵤ⁺⁻ - bᵤ⁺⁻
    Q⁺⁻ = Q + rᵥ⁺⁻ - bᵥ⁺⁻

    bᵤ⁻⁺ = zeros(Float64,T^2)
    bᵥ⁻⁺ = zeros(Float64,T^2)
    sᵤ⁻⁺ = Dᵤ⁻*Z - P + bᵤ⁻⁺
    sᵥ⁻⁺ = Dᵥ⁺*Z - Q + bᵥ⁻⁺
    rᵤ⁻⁺ = max(norm(sᵤ⁻⁺) - (4/α),0)*(sᵤ⁻⁺./(norm(sᵤ⁻⁺)+eps(Float64)))
    rᵥ⁻⁺ = max(norm(sᵥ⁻⁺) - (4/α),0)*(sᵥ⁻⁺./(norm(sᵥ⁻⁺)+eps(Float64)))
    bᵤ⁻⁺ = bᵤ⁻⁺ + Dᵤ⁻*Z - P - rᵤ⁻⁺
    bᵥ⁻⁺ = bᵥ⁻⁺ + Dᵥ⁺*Z - Q - rᵥ⁻⁺
    P⁻⁺ = P + rᵤ⁻⁺ - bᵤ⁻⁺
    Q⁻⁺ = Q + rᵥ⁻⁺ - bᵥ⁻⁺

    bᵤ⁻⁻ = zeros(Float64,T^2)
    bᵥ⁻⁻ = zeros(Float64,T^2)
    sᵤ⁻⁻ = Dᵤ⁻*Z - P + bᵤ⁻⁻
    sᵥ⁻⁻ = Dᵥ⁻*Z - Q + bᵥ⁻⁻
    rᵤ⁻⁻ = max(norm(sᵤ⁻⁻) - (4/α),0)*(sᵤ⁻⁻./(norm(sᵤ⁻⁻)+eps(Float64)))
    rᵥ⁻⁻ = max(norm(sᵥ⁻⁻) - (4/α),0)*(sᵥ⁻⁻./(norm(sᵥ⁻⁻)+eps(Float64)))
    bᵤ⁻⁻ = bᵤ⁻⁻ + Dᵤ⁻*Z - P - rᵤ⁻⁻
    bᵥ⁻⁻ = bᵥ⁻⁻ + Dᵥ⁻*Z - Q - rᵥ⁻⁻
    P⁻⁻ = P + rᵤ⁻⁻ - bᵤ⁻⁻
    Q⁻⁻ = Q + rᵥ⁻⁻ - bᵥ⁻⁻

    # Calculate bₜᵥ as per eq 60
    bₜᵥ = (transpose(Dᵤ⁺) * P⁺⁺) .+ (transpose(Dᵥ⁺) * Q⁺⁺) .+ (transpose(Dᵤ⁺) * P⁺⁻) .+ (transpose(Dᵥ⁻) * Q⁺⁻) .+ (transpose(Dᵤ⁻) * P⁻⁺) .+ (transpose(Dᵥ⁺) * Q⁻⁺) .+ (transpose(Dᵤ⁻) * P⁻⁻) .+ (transpose(Dᵥ⁻) * Q⁻⁻) .+ (Λ^2)*Z⁰
    for i = 1:max_iter
        cg!(Z, A, bₜᵥ, tol=10.0^-4, Pl = pl)

        # Recalculate b, s, r, P and Q in each combination of {+,-} as per equations 57, 61, 62
        sᵤ⁺⁺ = Dᵤ⁺*Z - P + bᵤ⁺⁺
        sᵥ⁺⁺ = Dᵥ⁺*Z - Q + bᵥ⁺⁺
        rᵤ⁺⁺ = max(norm(sᵤ⁺⁺) - (4/α),0)*(sᵤ⁺⁺./(norm(sᵤ⁺⁺)+eps(Float64)))
        rᵥ⁺⁺ = max(norm(sᵥ⁺⁺) - (4/α),0)*(sᵥ⁺⁺./(norm(sᵥ⁺⁺)+eps(Float64)))
        bᵤ⁺⁺ = bᵤ⁺⁺ + Dᵤ⁺*Z - P - rᵤ⁺⁺
        bᵥ⁺⁺ = bᵥ⁺⁺ + Dᵥ⁺*Z - Q - rᵥ⁺⁺
        P⁺⁺ = P + rᵤ⁺⁺ - bᵤ⁺⁺
        Q⁺⁺ = Q + rᵥ⁺⁺ - bᵥ⁺⁺

        sᵤ⁺⁻ = Dᵤ⁺*Z - P + bᵤ⁺⁻
        sᵥ⁺⁻ = Dᵥ⁻*Z - Q + bᵥ⁺⁻
        rᵤ⁺⁻ = max(norm(sᵤ⁺⁻) - (4/α),0)*(sᵤ⁺⁻./(norm(sᵤ⁺⁻)+eps(Float64)))
        rᵥ⁺⁻ = max(norm(sᵥ⁺⁻) - (4/α),0)*(sᵥ⁺⁻./(norm(sᵥ⁺⁻)+eps(Float64)))
        bᵤ⁺⁻ = bᵤ⁺⁻ + Dᵤ⁺*Z - P - rᵤ⁺⁻
        bᵥ⁺⁻ = bᵥ⁺⁻ + Dᵥ⁻*Z - Q - rᵥ⁺⁻
        P⁺⁻ = P + rᵤ⁺⁻ - bᵤ⁺⁻
        Q⁺⁻ = Q + rᵥ⁺⁻ - bᵥ⁺⁻

        sᵤ⁻⁺ = Dᵤ⁻*Z - P + bᵤ⁻⁺
        sᵥ⁻⁺ = Dᵥ⁺*Z - Q + bᵥ⁻⁺
        rᵤ⁻⁺ = max(norm(sᵤ⁻⁺) - (4/α),0)*(sᵤ⁻⁺./(norm(sᵤ⁻⁺)+eps(Float64)))
        rᵥ⁻⁺ = max(norm(sᵥ⁻⁺) - (4/α),0)*(sᵥ⁻⁺./(norm(sᵥ⁻⁺)+eps(Float64)))
        bᵤ⁻⁺ = bᵤ⁻⁺ + Dᵤ⁻*Z - P - rᵤ⁻⁺
        bᵥ⁻⁺ = bᵥ⁻⁺ + Dᵥ⁺*Z - Q - rᵥ⁻⁺
        P⁻⁺ = P + rᵤ⁻⁺ - bᵤ⁻⁺
        Q⁻⁺ = Q + rᵥ⁻⁺ - bᵥ⁻⁺

        sᵤ⁻⁻ = Dᵤ⁻*Z - P + bᵤ⁻⁻
        sᵥ⁻⁻ = Dᵥ⁻*Z - Q + bᵥ⁻⁻
        rᵤ⁻⁻ = max(norm(sᵤ⁻⁻) - (4/α),0)*(sᵤ⁻⁻./(norm(sᵤ⁻⁻)+eps(Float64)))
        rᵥ⁻⁻ = max(norm(sᵥ⁻⁻) - (4/α),0)*(sᵥ⁻⁻./(norm(sᵥ⁻⁻)+eps(Float64)))
        bᵤ⁻⁻ = bᵤ⁻⁻ + Dᵤ⁻*Z - P - rᵤ⁻⁻
        bᵥ⁻⁻ = bᵥ⁻⁻ + Dᵥ⁻*Z - Q - rᵥ⁻⁻
        P⁻⁻ = P + rᵤ⁻⁻ - bᵤ⁻⁻
        Q⁻⁻ = Q + rᵥ⁻⁻ - bᵥ⁻⁻

        # Recalculate bₜᵥ as per eq 60
        bₜᵥ = (transpose(Dᵤ⁺) * P⁺⁺) .+ (transpose(Dᵥ⁺) * Q⁺⁺) .+ (transpose(Dᵤ⁺) * P⁺⁻) .+ (transpose(Dᵥ⁻) * Q⁺⁻) .+ (transpose(Dᵤ⁻) * P⁻⁺) .+ (transpose(Dᵥ⁺) * Q⁻⁺) .+ (transpose(Dᵤ⁻) * P⁻⁻) .+ (transpose(Dᵥ⁻) * Q⁻⁻) .+ (Λ^2)*Z⁰
    end
    for j in eachindex(index)
        z[index[j]] = Z[j]
    end
    z = z.*mask
    z = rotl90(z)
    return z
end

@doc raw"""
```
NonConvex1(z::AbstractArray, β::Real = 0.5, λ::AbstractArray = fill(10.0^-6, size(z)), mask::AbstractArray = fill(1.0, size(z)), max_iter::Int = 100)
```
The first of two non-convex regularization methods proposed by Aujol, Durou and
Quéau. The same fidelity and regularization terms to minimize as [`Quadratic`](@ref)
are used, but the convexity of ``\Phi`` is sacrificed in order to gain better
behaviour around discontinuities and outliers.
# Output
`NonConvex1()` returns a NonConvex1 integrator which can then be called to run the
non-convex regularization method method on a gradient field.
# Details
The initial minimization problem is the same as in [`TotalVariation`](@ref)
as given below exept ``Phi`` is defined using the given function instead of the
``L_1`` norm.
```math
\begin{aligned}
&min\iint_{(u,v)\in \Omega}||\nabla z(u,v)-\mathbf{g}(u,v)||+\lambda(u,v)\left[z(u,v)-z^0(u,v)\right]^2dudv\\
&\text{where:}\\
&\Phi(s)=log(s^2+\beta^2)
\end{aligned}
```
The problem can then be discretized to form the following functional where ``D_{u,v}^{UV}``
is a ``2\times|\Omega|`` matrix formed from stacking the vectors of the finite
diferences.
```math
\begin{aligned}
E(\mathbf{z})&=\frac{1}{4}\sum_{(U,V)\in\{+,-\}^2}\sum_{(u,v)\in\Omega^{UV}}\Phi(||\nabla z_{u,v}-\mathbf{g}_{u,v}||)+\sum_{(u,v)\in\Omega}\lambda_{u,v}\left[z_{u,v}-z^0_{u,v}\right]^2\\
&=\frac{1}{4}\sum_{(U,V)\in\{+,-\}^2}\sum_{(u,v)\in\Omega^{UV}}\Phi(||D_{u,v}^{UV}\mathbf{z}-\mathbf{g}_{u,v}||)+||\Lambda\left(\mathbf{z}-\mathbf{z}^0\right)||^2\\
&=f(\mathbf{z})+g(\mathbf{z})
\end{aligned}
```
As ``f(\mathbf{z})`` is smooth but not convex and ``g(\mathbf{z})`` is convex
the iPiano algorithm is then used to iterativly solve the minimization of the
functional such that:
```math
\mathbf{z}^{(k+1)}=(I+\alpha_1\partial g)^{-1}(\mathbf{z}^{(k)}-\alpha_1\nabla f(\mathbf{z}^{(k)})+\alpha_2(\mathbf{z}^{(k)}-z^{(k+1)}))
```
where ``(I+\alpha_1\partial g)^{-1}`` is a proximal operator defined as:
```math
(I+\alpha_1\partial g)^{-1}(\mathbf{\hat{x}})=(I+2\alpha_1\Lambda^2)^{-1}(\mathbf{\hat{x}}+2\alpha_1\Lambda\mathbf{z}^0)
```
Using the definition of ``\Phi`` given above, the derivative of ``f(\mathbf{z})``
can be computed as below:
```math
\nabla f(\mathbf{z})=\frac{1}{4}\sum_{(U,V)\in\{+,-\}^2}\sum_{(u,v)\in\Omega^{UV}}\frac{D_{u,v}^{UV^\top}(D_{u,v}^{UV}\mathbf{z}-\mathbf{g}_{u,v})}{||D_{u,v}^{UV}\mathbf{z}-\mathbf{g}_{u,v}||^2+\beta^2}
```
The values of ``\alpha_1`` and ``\aplha_2`` control the step size of the algorith
and are chosen such that ``\alpha_2`` is fixed at 0.8 while ``\alpha_1`` is chosen
using the lazy backtracking method, which uses a Lipschitz constant to caculate
a suitably small step size at each step using the below relationship where ``\eta>1``
is a constant:
```math
\begin{aligned}
&\alpha_1 < 2(1-\alpha_2)/L_n\\
&\text{where:}\\
&L_k\in\{L_{k-1},\eta L_{k-1},\eta^2L_{k-1},...\}\\
&\text{such that it is minimal and satisfies:}\\
&f(x^{(k+1)})\le f(x^{(k)}) + \left\langle\nabla f(x^{(k)}), x^{(k+1)}-x^{(k)}\right\rangle+\frac{L_k}{2}||x^{(k+1)}-x^{(k)}||^2
\end{aligned}
```
This method, being non-convex, highly relies on a good initial solution and will
often only provide a minimal improvment to the solution. A bad initial solution
will produce a final solution which does not resemble the surface under most
conditions.
# Parameters
## `z`:
An `AbstractArray` which defines the value of ``z^0`` the initial solution and
prior to be used in the regulization term.
## `β`:
A `Real` which acts as a hyper-parameter to the function. Large values for β will
produce smoother functions but loose discontinuities. Smaller value will preserve
discontinuities but lead to staircassing in the solution. Defults to `0.5`.
## `λ`:
An `AbstractArray` the same size as `z` defulting to ``10.0^-6`` everywhere.
This defines theregulization weight at each point. Large valueas will force the
algorithm to keep the solution near to ``z^0`` at that position. Can be used to
keep the solution near the initial solution or guide the solution to a certian
known value at points (i.e. known maxima and minima). This value should be set
uniformly small otherwise.
## `mask`:
An `AbstractArray` the same size as `z`, which guides the algorithm as to where
the valid domain is. Values of `1` will be in the domain ``\Omega`` while other
values will be ignored and set to ``z^0``. This can be used to
integrate over sub-domain or to segment the domain into parts. The `gen_mask()`
funtion can be used to generate a mask which will remove non-integrable regions
dramatically improving the solution under most condition at the cost of not
integrating the entire solution.
# Example
The following example demonstrates the use of the `NonConvex1` integrator.
```julia
using ShapeFromShading, Makie

# Generate synthetic gradients
p, q = synthetic_gradient(Prism(75), img_size = 151)

# Create a NonConvex1() integrator
nonConvex1 = NonConvex1(z=zeros(size(p)), β=0.5)
nonConvex1Init = NonConvex1(z=Horn()(p,q), β=0.5)

# Calculate the heightmap from the gradients
Z = nonConvex1(p, q)
Z2 = nonConvex1Init(p, q)

# Normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(Z)
Z2 = Z2./maximum(Z2)

# Display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
vbox(surface(r, r, Z), surface(r, r, Z2))
```
# Reference
[1] Y. Quéau, J. Durou and J. Aujol, "Variational Methods for Normal Integration", Journal of Mathematical Imaging and Vision, vol. 60, no. 4, pp. 609-632, 2017. [doi: 10.1007/s10851-017-0777-6](https://doi.org/10.1007/s10851-017-0777-6 )
"""
function (scheme::NonConvex1)(pIn::AbstractArray, qIn::AbstractArray)
    z⁰ = copy(scheme.z)
    λ = copy(scheme.λ)
    mask = copy(scheme.mask)
    β = copy(scheme.β)
    max_iter = copy(scheme.max_iter)
    z = copy(z⁰)
    mask = rotr90(mask)
    λ = rotr90(λ)
    p = rotr90(copy(pIn)).*mask
    q = rotr90(copy(qIn)).*mask
    T = first(size(z))
    index = CartesianIndices(z)
    index = reshape(index, length(index),1)

    Dᵤ⁺ = spzeros(Float64,T^2,T^2)
    Dᵤ⁻ = spzeros(Float64,T^2,T^2)
    Dᵥ⁺ = spzeros(Float64,T^2,T^2)
    Dᵥ⁻ = spzeros(Float64,T^2,T^2)
    gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻, T, mask)
    Λ = zeros(Float64,T^2,T^2)
    P = zeros(Float64,T^2)
    Q = zeros(Float64,T^2)
    Z⁰ = zeros(Float64,T^2)
    Z = zeros(Float64,T^2)

    # Vectorize inputs
    for j in CartesianIndices(z)
        Λ[j[1]+(j[2]-1)*T,j[1]+(j[2]-1)*T] = sqrt(λ[j])
        P[j[1]+(j[2]-1)*T] = p[j]
        Q[j[1]+(j[2]-1)*T] = q[j]
        Z⁰[j[1]+(j[2]-1)*T] = z⁰[j]
        Z[j[1]+(j[2]-1)*T] = z[j]
    end
    Λ = Diagonal(Λ)
    Zold = copy(Z)
    Ztemp = copy(Z)
    α₁ = 0.8
    α₂ = 0.8

    # Initialize Lₙ (Lipschitz constant) such that lazy backtracking starts with α₁ = α₂.
    # This is mostly arbitary and value will approach needed value over time.
    Lₙ = (18.0/α₁)*(1-α₂)
    η = 1.1
    i = 1
    while i <= max_iter && (norm(Z - Zold) > 10.0^-4.0 || i == 1)
        # Calculate f₁, ∇f₁ and previous f₁ (f₁ⁿ⁻¹ )
        ∇f₁ = ((transpose(Dᵤ⁺)*(Dᵤ⁺*Z-P))/(norm(Dᵤ⁺*Z-P)^2+β^2)) + ((transpose(Dᵥ⁺)*(Dᵥ⁺*Z-Q))/(norm(Dᵥ⁺*Z-Q)^2+β^2)) + ((transpose(Dᵤ⁻)*(Dᵤ⁻*Z-P))/(norm(Dᵤ⁻*Z-P)^2+β^2)) + ((transpose(Dᵥ⁻)*(Dᵥ⁻*Z-Q))/(norm(Dᵥ⁻*Z-Q)^2+β^2))
        f₁ⁿ⁻¹ = log(norm(Dᵤ⁺*Zold-P)^2+β^2) + log(norm(Dᵥ⁺*Zold-Q)^2+β^2) + log(norm(Dᵤ⁻*Zold-P)^2+β^2) + log(norm(Dᵥ⁻*Zold-Q)^2+β^2)
        f₁ = log(norm(Dᵤ⁺*Z-P)^2+β^2) + log(norm(Dᵥ⁺*Z-Q)^2+β^2) + log(norm(Dᵤ⁻*Z-P)^2+β^2) + log(norm(Dᵥ⁻*Z-Q)^2+β^2)

        # Calculate lazy backtracking condition as per ipiano algortihm 4
        while f₁ > f₁ⁿ⁻¹  + dot(∇f₁, Z - Zold) + (Lₙ/2)*norm(Z - Zold)
            Lₙ = Lₙ*η
        end

        # Calculate new α₁ as per ipiano algortihm 4
        α₁ = (2*(1-α₂)/Lₙ)*0.9

        # Calcualte x̂ as per eq 69 then calculate new Z as per eq 71
        x̂ = Z - α₁.*∇f₁ + α₂.*(Z-Zold)
        Ztemp = inv(I+(2*α₁).*Λ^2)*(x̂+(2*α₁).*Λ*Z⁰)
        copyto!(Zold, Z)
        copyto!(Z, Ztemp)
        i += 1
    end
    println("Converged after ", i-1, " iterations")
    for j in eachindex(index)
        z[index[j]] = Z[j]
    end
    z = rotl90(z)
    return z
end

@doc raw"""
```
NonConvex2(z::AbstractArray, γ::Real = 1.0, λ::AbstractArray = fill(10.0^-6, size(z, mask::AbstractArray = fill(1.0, size(z)), max_iter::Int = 100)
```
The second of two non-convex regularization methods proposed by Aujol, Durou and
Quéau. The same as [`NonConvex1`](@ref) exept ``\Phi`` has been replaced with a
different non-convex function.
# Output
`NonConvex2()` returns a NonConvex2 integrator which can then be called to run the
non-convex regularization method method on a gradient field.
# Details
The implimentation is the same as the one in [`TotalVariation`](@ref)
as given below exept ``Phi`` is defined as below:
```math
\Phi(s)=\frac{s^2}{s^2+\gamma^2}
```
Using the definition of ``\Phi`` given above, the derivative of ``f(\mathbf{z})``
can be computed as below:
```math
\nabla f(\mathbf{z})=\frac{1}{4}\sum_{(U,V)\in\{+,-\}^2}\sum_{(u,v)\in\Omega^{UV}}\frac{\gamma^2D_{u,v}^{UV^\top}(D_{u,v}^{UV}\mathbf{z}-\mathbf{g}_{u,v})}{()||D_{u,v}^{UV}\mathbf{z}-\mathbf{g}_{u,v}||^2+\gamma^2)^2}
```
This method, being non-convex, highly relies on a good initial solution and will
often only provide a minimal improvment to the solution. A bad initial solution
will produce a final solution which does not resemble the surface under most
conditions.
# Parameters
## `z`:
An `AbstractArray` which defines the value of ``z^0`` the initial solution and
prior to be used in the regulization term.
## `γ`:
A `Real` which acts as a hyper-parameter to the function. Large values for γ will
produce smoother functions but loose discontinuities. Smaller value will preserve
discontinuities but lead to staircassing in the solution. Defults to `1.0`.
## `λ`:
An `AbstractArray` the same size as `z` defulting to ``10.0^-6`` everywhere.
This defines theregulization weight at each point. Large valueas will force the
algorithm to keep the solution near to ``z^0`` at that position. Can be used to
keep the solution near the initial solution or guide the solution to a certian
known value at points (i.e. known maxima and minima). This value should be set
uniformly small otherwise.
## `mask`:
An `AbstractArray` the same size as `z`, which guides the algorithm as to where
the valid domain is. Values of `1` will be in the domain ``\Omega`` while other
values will be ignored and set to ``z^0``. This can be used to
integrate over sub-domain or to segment the domain into parts. The `gen_mask()`
funtion can be used to generate a mask which will remove non-integrable regions
dramatically improving the solution under most condition at the cost of not
integrating the entire solution.
# Example
The following example demonstrates the use of the `NonConvex2` integrator.
```julia
using ShapeFromShading, Makie

# Generate synthetic gradients
p, q = synthetic_gradient(Prism(75), img_size = 151)

# Create a NonConvex2() integrator
nonConvex2 = NonConvex2(z=zeros(size(p)), γ=1.0)
nonConvex2Init = NonConvex2(z=Horn()(p,q), γ=1.0)

# Calculate the heightmap from the gradients
Z = nonConvex2(p, q)
Z2 = nonConvex2Init(p, q)

# Normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(Z)
Z2 = Z2./maximum(Z2)

# Display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
vbox(surface(r, r, Z), surface(r, r, Z2))
```
# Reference
[1] Y. Quéau, J. Durou and J. Aujol, "Variational Methods for Normal Integration", Journal of Mathematical Imaging and Vision, vol. 60, no. 4, pp. 609-632, 2017. [doi: 10.1007/s10851-017-0777-6](https://doi.org/10.1007/s10851-017-0777-6 )
"""
function (scheme::NonConvex2)(pIn::AbstractArray, qIn::AbstractArray)
    z⁰ = copy(scheme.z)
    λ = copy(scheme.λ)
    mask = copy(scheme.mask)
    γ = copy(scheme.γ)
    max_iter = copy(scheme.max_iter)
    z = copy(z⁰)
    mask = rotr90(mask)
    λ = rotr90(λ)
    p = rotr90(copy(pIn)).*mask
    q = rotr90(copy(qIn)).*mask
    T = first(size(z))
    index = CartesianIndices(z)
    index = reshape(index, length(index),1)

    Dᵤ⁺ = spzeros(Float64,T^2,T^2)
    Dᵤ⁻ = spzeros(Float64,T^2,T^2)
    Dᵥ⁺ = spzeros(Float64,T^2,T^2)
    Dᵥ⁻ = spzeros(Float64,T^2,T^2)
    gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻, T, mask)
    Λ = zeros(Float64,T^2,T^2)
    P = zeros(Float64,T^2)
    Q = zeros(Float64,T^2)
    Z⁰ = zeros(Float64,T^2)
    Z = zeros(Float64,T^2)

    # Vectorize inputs
    for j in CartesianIndices(z)
        Λ[j[1]+(j[2]-1)*T,j[1]+(j[2]-1)*T] = sqrt(λ[j])
        P[j[1]+(j[2]-1)*T] = p[j]
        Q[j[1]+(j[2]-1)*T] = q[j]
        Z⁰[j[1]+(j[2]-1)*T] = z⁰[j]
        Z[j[1]+(j[2]-1)*T] = z[j]
    end
    Λ = Diagonal(Λ)
    Zold = copy(Z)
    Ztemp = copy(Z)
    α₁ = 0.8
    α₂ = 0.8

    # Initialize Lₙ (Lipschitz constant) such that lazy backtracking starts with α₁ = α₂.
    # This is mostly arbitary and value will approach needed value over time.
    Lₙ = (18.0/α₁)*(1-α₂)
    η = 1.1
    i = 1
    while i <= max_iter && (norm(Z - Zold) > 10.0^-4.0 || i == 1)
        # Calculate f₁, ∇f₁ and previous f₁ (f₁ⁿ⁻¹ )
        ∇f₁ = (((γ^2)*transpose(Dᵤ⁺)*(Dᵤ⁺*Z-P))/(norm(Dᵤ⁺*Z-P)^2+γ^2)^2) + (((γ^2)*transpose(Dᵥ⁺)*(Dᵥ⁺*Z-Q))/(norm(Dᵥ⁺*Z-Q)^2+γ^2)^2) + (((γ^2)*transpose(Dᵤ⁻)*(Dᵤ⁻*Z-P))/(norm(Dᵤ⁻*Z-P)^2+γ^2)^2) + (((γ^2)*transpose(Dᵥ⁻)*(Dᵥ⁻*Z-Q))/(norm(Dᵥ⁻*Z-Q)^2+γ^2)^2)
        f₁ⁿ⁻¹ = (norm(Dᵤ⁺*Zold-P)^2)/(norm(Dᵤ⁺*Zold-P)^2+γ^2) + (norm(Dᵥ⁺*Zold-Q)^2)/(norm(Dᵥ⁺*Zold-Q)^2+γ^2) + (norm(Dᵤ⁻*Zold-P)^2)/(norm(Dᵤ⁻*Zold-P)^2+γ^2) + (norm(Dᵥ⁻*Zold-Q)^2)/(norm(Dᵥ⁻*Zold-Q)^2+γ^2)
        f₁ = (norm(Dᵤ⁺*Z-P)^2)/(norm(Dᵤ⁺*Z-P)^2+γ^2) + (norm(Dᵥ⁺*Z-Q)^2)/(norm(Dᵥ⁺*Z-Q)^2+γ^2) + (norm(Dᵤ⁻*Z-P)^2)/(norm(Dᵤ⁻*Z-P)^2+γ^2) + (norm(Dᵥ⁻*Z-Q)^2)/(norm(Dᵥ⁻*Z-Q)^2+γ^2)

        # Calculate lazy backtracking condition as per ipiano algortihm 4
        while f₁ > f₁ⁿ⁻¹ + dot(∇f₁, Z - Zold) + (Lₙ/2)*norm(Z - Zold)
            Lₙ = Lₙ*η
        end

        # Calculate new α₁ as per ipiano algortihm 4
        α₁ = (2*(1-α₂)/Lₙ)*0.9

        # Calcualte x̂ as per eq 69 then calculate new Z as per eq 71
        x̂ = Z - α₁.*∇f₁ + α₂.*(Z-Zold)
        Ztemp = inv(I+(2*α₁).*Λ^2)*(x̂+(2*α₁).*Λ*Z⁰)
        copyto!(Zold, Z)
        copyto!(Z, Ztemp)
        i += 1
    end
    println("Converged after ", i-1, " iterations")
    for j in eachindex(index)
        z[index[j]] = Z[j]
    end
    z = rotl90(z)
    return z
end

@doc raw"""
```
AnisotropicDiffusion(z::AbstractArray, μ::Real = 5.0, ν::Real = 10., λ::AbstractArray = fill(10.0^-6, size(z)), mask::AbstractArray = fill(1.0, size(z)), max_iter::Int = 10)
```
Defines the anisotropic diffusion method as proposed by Aujol, Durou and Quéau.
It utilizes an anisotropic diffusion like process to create a weighted version
of the least squares problem solved in the [`Quadratic`](@ref) method.
# Output
`AnisotropicDiffusion()` returns a AnisotropicDiffusion integrator which can then be called to run the
anisotropic diffusion method method on a gradient field.
# Details
The minimization problem from [`Quadratic`](@ref) is modified with the addition
weighting term ``W(u,v)`` to reach the following minimization problem.
```math
min\iint_{(u,v)\in \Omega}||W(u,v)[\nabla z(u,v)-\mathbf{g}(u,v)]||^2+\lambda(u,v)\left[z(u,v)-z^0(u,v)\right]^2dudv
```
The weighting term can be defined as below where ``\mu`` and ``\nu`` are parameters
which control the impact of gradient of the reconstruction at each iteration and
the input gradients.
```math
W(u,v)=\frac{1}{\sqrt{\left(\frac{||\nabla z(u,v)||}{\mu}\right)^2+1}}\begin{bmatrix}\frac{1}{\sqrt{1+\left(\frac{p(u,v)}{\nu}\right)^2}} &0\\ 0 &\frac{1}{\sqrt{1+\left(\frac{q(u,v)}{\nu}\right)^2}}\end{bmatrix}
```
This lends itself to the creation of the terms ``A^{UV}`` and ``B^{UV}`` which
are ``|\Omega|\times|\Omega|`` diagonal matrices, respectively,  contianing the
values of;
```math
    \frac{1}{\sqrt{1+\left(\frac{p(u,v)}{\nu}\right)^2}\sqrt{\frac{(\partial^U_uz_{u,v})^2+(\partial^V_vz_{u,v})^2}{\mu^2}+1}}
```
and
```math
    \frac{1}{\sqrt{1+\left(\frac{q(u,v)}{\nu}\right)^2}\sqrt{\frac{(\partial^U_uz_{u,v})^2+(\partial^V_vz_{u,v})^2}{\mu^2}+1}}
```
where ``(U,V)\in\{+,-\}^2``. Using these definititions the original minimization
problem can be rewriten as the following iterative scheme:
```math
\begin{aligned}
z^{(k+1)}=&min\frac{1}{4}\sum_{(U,V)\in\{+,-\}^2}\left\{||A^{UV}(D^U_u\mathbf{z}-\mathbf{p})||^2+||B^{UV}(D^V_v\mathbf{z}-\mathbf{q})||^2\right\}\\
&+||\Lambda(\mathbf{z}-\mathbf{z}^0)||
\end{aligned}
```
which is then solved using Cholesky factorization.
# Parameters
## `z`:
An `AbstractArray` which defines the value of ``z^0`` the initial solution and
prior to be used in the regulization term.
## `μ`:
A `Real` which acts as a hyper-parameter to the function. Allows for the tuning
of discontinuities. Defults to `5.0`.
## `ν`:
A `Real` which acts as a hyper-parameter to the function. Allows for the tuning
of discontinuities. Unlike `μ` is has a minimal impact on the solution and generally
does not need adjusting unless an extreme value for `μ` is used where it can help
to balance out the two terms. Defults to `10.0`.
## `λ`:
An `AbstractArray` the same size as `z` defulting to ``10.0^-6`` everywhere.
This defines theregulization weight at each point. Large valueas will force the
algorithm to keep the solution near to ``z^0`` at that position. Can be used to
keep the solution near the initial solution or guide the solution to a certian
known value at points (i.e. known maxima and minima). This value should be set
uniformly small otherwise.
## `mask`:
An `AbstractArray` the same size as `z`, which guides the algorithm as to where
the valid domain is. Values of `1` will be in the domain ``\Omega`` while other
values will be ignored and set to ``z^0``. This can be used to
integrate over sub-domain or to segment the domain into parts. The `gen_mask()`
funtion can be used to generate a mask which will remove non-integrable regions
dramatically improving the solution under most condition at the cost of not
integrating the entire solution.
# Example
The following example demonstrates the use of the `AnisotropicDiffusion` integrator.
```julia
using ShapeFromShading, Makie

# Generate synthetic gradients
p, q = synthetic_gradient(Prism(75), img_size = 151)

# Create a AnisotropicDiffusion() integrator
anisotropicDiffusion = AnisotropicDiffusion(z=zeros(size(p)), γ=1.0)
anisotropicDiffusionInit = AnisotropicDiffusion(z=Horn()(p,q), γ=1.0)

# Calculate the heightmap from the gradients
Z = anisotropicDiffusion(p, q)
Z2 = anisotropicDiffusionInit(p, q)

# Normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(Z)
Z2 = Z2./maximum(Z2)

# Display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
vbox(surface(r, r, Z), surface(r, r, Z2))
```
# Reference
[1] Y. Quéau, J. Durou and J. Aujol, "Variational Methods for Normal Integration", Journal of Mathematical Imaging and Vision, vol. 60, no. 4, pp. 609-632, 2017. [doi: 10.1007/s10851-017-0777-6](https://doi.org/10.1007/s10851-017-0777-6 )
"""
function (scheme::AnisotropicDiffusion)(pIn::AbstractArray, qIn::AbstractArray)
    λ = copy(scheme.λ)
    z⁰ = copy(scheme.z)
    mask = copy(scheme.mask)
    ν = copy(scheme.ν)
    μ = copy(scheme.μ)
    max_iter = copy(scheme.max_iter)
    z = copy(z⁰)
    mask = rotr90(mask)
    λ = rotr90(λ)
    p = rotr90(copy(pIn)).*mask
    q = rotr90(copy(qIn)).*mask
    T = first(size(z))
    index = CartesianIndices(z)
    index = reshape(index, length(index),1)

    Dᵤ⁺ = spzeros(Float64,T^2,T^2)
    Dᵤ⁻ = spzeros(Float64,T^2,T^2)
    Dᵥ⁺ = spzeros(Float64,T^2,T^2)
    Dᵥ⁻ = spzeros(Float64,T^2,T^2)
    gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻,T,mask)
    Λ = zeros(Float64,T^2,T^2)
    P = zeros(Float64,T^2)
    Q = zeros(Float64,T^2)
    Z⁰ = zeros(Float64,T^2)
    Z = zeros(Float64,T^2)

    # Vectorize inputs
    for j in CartesianIndices(z)
        Λ[j[1]+(j[2]-1)*T,j[1]+(j[2]-1)*T] = sqrt(λ[j])
        P[j[1]+(j[2]-1)*T] = p[j]
        Q[j[1]+(j[2]-1)*T] = q[j]
        Z⁰[j[1]+(j[2]-1)*T] = z⁰[j]
        Z[j[1]+(j[2]-1)*T] = z[j]
    end
    Λ = sparse(Diagonal(Λ))

    # Calcualte A and B for each conination in {+,-} as per eq 90 and 91
    A⁺⁺ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
    A⁻⁺ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
    A⁺⁻ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))
    A⁻⁻ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))
    B⁺⁺ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
    B⁻⁺ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
    B⁺⁻ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))
    B⁻⁻ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))

    # Calculate A and b as derived from eq 92 in form Az=b
    A = (A⁺⁺*Dᵤ⁺)*transpose(A⁺⁺*Dᵤ⁺) + (B⁺⁺*Dᵥ⁺)*transpose(B⁺⁺*Dᵥ⁺) + (A⁺⁻*Dᵤ⁺)*transpose(A⁺⁻*Dᵤ⁺) + (B⁺⁻*Dᵥ⁻)*transpose(B⁺⁻*Dᵥ⁻) + (A⁻⁺*Dᵤ⁻)*transpose( A⁻⁺*Dᵤ⁻) + (B⁻⁺*Dᵥ⁺)*transpose(B⁻⁺*Dᵥ⁺) + (A⁻⁻*Dᵤ⁻)*transpose(A⁻⁻*Dᵤ⁻) + (B⁻⁻*Dᵥ⁻)*transpose(B⁻⁻*Dᵥ⁻)
    A = A/4.0 + Λ^2
    b = (transpose(A⁺⁺*Dᵤ⁺)*A⁺⁺)*P + (transpose(B⁺⁺*Dᵥ⁺)*B⁺⁺)*Q + (transpose(A⁺⁻*Dᵤ⁺)*A⁺⁻)*P + (transpose(B⁺⁻*Dᵥ⁻)*B⁺⁻)*Q + (transpose(A⁻⁺*Dᵤ⁻)*A⁻⁺)*P + (transpose(B⁻⁺*Dᵥ⁺)*B⁻⁺)*Q + (transpose(A⁻⁻*Dᵤ⁻)*A⁻⁻)*P + (transpose(B⁻⁻*Dᵥ⁻)*B⁻⁻)*Q
    b = b/4.0 + (Λ^2)*Z⁰
    A₀ = factorize(A)
    for i = 1:max_iter
        # Solve system at iteration k
        Z = A₀\b

        # Recalcualte A and B for each conination in {+,-} as per eq 90 and 91
        A⁺⁺ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
        A⁻⁺ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
        A⁺⁻ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))
        A⁻⁻ = Diagonal(1.0./(sqrt.(1.0.+(P./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))
        B⁺⁺ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
        B⁻⁺ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁺*Z).^2)./(μ^2) .+ 1)))
        B⁺⁻ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁺*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))
        B⁻⁻ = Diagonal(1.0./(sqrt.(1.0.+(Q./ν).^2) .* sqrt.(((Dᵤ⁻*Z).^2 .+ (Dᵥ⁻*Z).^2)./(μ^2) .+ 1)))

        # Recalculate A and b as derived from eq 92 in form Az=b
        A = (A⁺⁺*Dᵤ⁺)*transpose(A⁺⁺*Dᵤ⁺) + (B⁺⁺*Dᵥ⁺)*transpose(B⁺⁺*Dᵥ⁺) + (A⁺⁻*Dᵤ⁺)*transpose(A⁺⁻*Dᵤ⁺) + (B⁺⁻*Dᵥ⁻)*transpose(B⁺⁻*Dᵥ⁻) + (A⁻⁺*Dᵤ⁻)*transpose( A⁻⁺*Dᵤ⁻) + (B⁻⁺*Dᵥ⁺)*transpose(B⁻⁺*Dᵥ⁺) + (A⁻⁻*Dᵤ⁻)*transpose(A⁻⁻*Dᵤ⁻) + (B⁻⁻*Dᵥ⁻)*transpose(B⁻⁻*Dᵥ⁻)
        A = A/4.0 + Λ^2
        b = (transpose(A⁺⁺*Dᵤ⁺)*A⁺⁺)*P + (transpose(B⁺⁺*Dᵥ⁺)*B⁺⁺)*Q + (transpose(A⁺⁻*Dᵤ⁺)*A⁺⁻)*P + (transpose(B⁺⁻*Dᵥ⁻)*B⁺⁻)*Q + (transpose(A⁻⁺*Dᵤ⁻)*A⁻⁺)*P + (transpose(B⁻⁺*Dᵥ⁺)*B⁻⁺)*Q + (transpose(A⁻⁻*Dᵤ⁻)*A⁻⁻)*P + (transpose(B⁻⁻*Dᵥ⁻)*B⁻⁻)*Q
        b = b/4.0 + (Λ^2)*Z⁰
        A₀ = factorize(A)
    end

    for j in eachindex(index)
        z[index[j]] = Z[j]
    end
    z = z.*mask
    z = rotl90(z)
    return z
end

@doc raw"""
```
MumfordShah(z::AbstractArray, μ::Real = 10.0, ϵ::Real = 0.1, λ::AbstractArray = fill(10.0^-6, size(z, mask::AbstractArray = fill(1.0, size(z)), max_iter::Int = 50)
```
Defines the Mumford-Shah method, the final method proposed by Aujol, Durou and
Quéau. It utilizes an adepted version of Mumford and Shah functional to provide
a good reconstruction in the presence of discontinuities at the const of a longer
runtime then other methods.
!!! warning
    This function can take several minutes to run on grids larger then 64X64.
# Output
`MumfordShah()` returns a MumfordShah integrator which can then be called to run the
Mumford-Shah method method on a gradient field.
# Details
This method involves modifying the minimization problem given by the Mumford and
Shah functional given below:
```math
min\;\mu\iint_{(u,v)\in\Omega\backslash K}||\nabla z(u,v)||dudv+\int_Kd\sigma\\+\lambda\iint_{(u,v)\in\Omega\backslash K}[z(u,v)-z^0(u,v)]^2dudv
```
Where ``K`` is the set of discontinuities. This can then be adapted to our problem
by utilizing a Ambrosio-Tortelli approximation to achaive the following functional
where ``\epsilon\to0``.
```math
\begin{aligned}
E(z)=\mu&\iint_{(u,v)\in\Omega}w(u,v)^2||\nabla z(u,v)-\mathbf{g}(u,v)||^2dudv\\
+&\iint_{(u,v)}\left[\epsilon||\nabla w(u,v)||^2+\frac{1}{4\epsilon}(w(u,v)-1)^2\right]dudv\\
+&\iint_{(u,v)}\lambda(u,v)\left[z(u,v)-z^0(u,v)\right]dudv
\end{aligned}
```
This leads to the following discretization of the functional where ``\mathbf{w}^{+/-}_{u/v}``
are thevecotr of weights in each forward and backward direction and the matrix
``W^{+/-}_{u/v}`` is the diagonal matrix formed from this vector.
```math
\begin{aligned}
E(\mathbf{z},\mathbf{w}^{+}_{u},\mathbf{w}^{-}_{u},\mathbf{w}^{+}_{v},\mathbf{w}^{-}_{v})&=\frac{\mu}{2}\Big(||W^+_u(D^+_u\mathbf{z}-\mathbf{p})||^2+||W^-_u(D^-_u\mathbf{z}-\mathbf{p})||^2\\
&+||W^+_v(D^+_v\mathbf{z}-\mathbf{q})||^2+||W^-_v(D^-_v\mathbf{z}-\mathbf{q})||^2\Big)\\
&+\frac{\epsilon}{2}\Big(||D^+_u\mathbf{w}^{+}_{u}||^2+||D^-_u\mathbf{w}^{-}_{u}||^2+||D^+_v\mathbf{w}^{+}_{v}||^2\\
&+||D^-_v\mathbf{w}^{-}_{v}||^2\Big)+\frac{1}{8\epsilon}\Big(||\mathbf{w}^{+}_{u}-\mathbf{1}||^2+||\mathbf{w}^{-}_{u}-\mathbf{1}||^2\\
&+||\mathbf{w}^{+}_{v}-\mathbf{1}||^2+||\mathbf{w}^{-}_{v}-\mathbf{1}||^2\Big)+||\Lambda(\mathbf{z}-\mathbf{z}^0)||^2
\end{aligned}
```
This is then solved with a conjugate gradient algorithm and an alternating optimization
scheme at each step where the updates are found using the relationships below:
```math
\begin{gathered}
\mathbf{z}^{(k+1)}=min\;E(\mathbf{z}^{(k)},\mathbf{w}^{+(k)}_{u},\mathbf{w}^{-(k)}_{u},\mathbf{w}^{+(k)}_{v},\mathbf{w}^{-(k)}_{v})\\
\mathbf{w}^{+/-(k+1)}_{u/v}=min\;E(\mathbf{z}^{(k+1},\mathbf{w}^{+(k)}_{u},\mathbf{w}^{-(k)}_{u},\mathbf{w}^{+(k)}_{v},\mathbf{w}^{-(k)}_{v})\\
\end{gathered}
```
This method can produce good results if the parameters are appropriotly tuned.
If the parameters are too large the solution will suffer heavily from staircasing
artifacts while setting it too small will result in a smooth solution. Even if the
value is chosen correctly the algorithm tends to overfit the final solution to the
discontinuities and they will extend into parts of the solution which are actually
smooth.
# Parameters
## `z`:
An `AbstractArray` which defines the value of ``z^0`` the initial solution and
prior to be used in the regulization term.
## `ϵ`:
A `Real` which acts as a hyper-parameter to the function. Controls how the final
solution will converge. Large values will lead to staircasing while small values
will over-smooth the surface. Must be relativly small to achieve convergence to a
solution. Defults to `0.1`.
## `μ`:
A `Real` which acts as a hyper-parameter to the function. This value controls the
smoothness of the final solution. Large values will lead to staircasing while
small values will lead to over-smoothed solutions. Defults to `10.0`.
## `λ`:
An `AbstractArray` the same size as `z` defulting to ``10.0^-6`` everywhere.
This defines the regulization weight at each point. Large values will force the
algorithm to keep the solution near to ``z^0`` at that position. Can be used to
keep the solution near the initial solution or guide the solution to a certian
known value at points (i.e. known maxima and minima). This value should be set
uniformly small otherwise.
## `mask`:
An `AbstractArray` the same size as `z`, which guides the algorithm as to where
the valid domain is. Values of `1` will be in the domain ``\Omega`` while other
values will be ignored and set to ``z^0``. This can be used to
integrate over sub-domain or to segment the domain into parts. The `gen_mask()`
funtion can be used to generate a mask which will remove non-integrable regions
dramatically improving the solution under most condition at the cost of not
integrating the entire solution.
# Example
The following example demonstrates the use of the `MumfordShah` integrator.
```julia
using ShapeFromShading, Makie

# Generate synthetic gradients
p, q = synthetic_gradient(Prism(75), img_size = 151)

# Create a MumfordShah() integrator
mumfordShah = MumfordShah(z=zeros(size(p)), γ=1.0)
mumfordShahInit = MumfordShah(z=Horn()(p,q), γ=1.0)

# Calculate the heightmap from the gradients
Z = mumfordShah(p, q)
Z2 = mumfordShahInit(p, q)

# Normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(Z)
Z2 = Z2./maximum(Z2)

# Display using Makie (Note: Makie can often take several minutes first time)
r = 0.0:0.1:4
vbox(surface(r, r, Z), surface(r, r, Z2))
```
# Reference
[1] Y. Quéau, J. Durou and J. Aujol, "Variational Methods for Normal Integration", Journal of Mathematical Imaging and Vision, vol. 60, no. 4, pp. 609-632, 2017. [doi: 10.1007/s10851-017-0777-6](https://doi.org/10.1007/s10851-017-0777-6 )
"""
function (scheme::MumfordShah)(pIn::AbstractArray, qIn::AbstractArray)
    z⁰ = copy(scheme.z)
    λ = copy(scheme.λ)
    mask = copy(scheme.mask)
    ϵ = copy(scheme.ϵ)
    μ = copy(scheme.μ)
    max_iter = copy(scheme.max_iter)
    z = copy(z⁰)
    mask = rotr90(mask)
    λ = rotr90(copy(λ))
    p = rotr90(copy(pIn)).*mask
    q = rotr90(copy(qIn)).*mask
    T = first(size(z))
    index = CartesianIndices(z)
    index = reshape(index, length(index),1)

    Dᵤ⁺ = spzeros(Float64,T^2,T^2)
    Dᵤ⁻ = spzeros(Float64,T^2,T^2)
    Dᵥ⁺ = spzeros(Float64,T^2,T^2)
    Dᵥ⁻ = spzeros(Float64,T^2,T^2)
    gen_matrix(Dᵤ⁺, Dᵤ⁻, Dᵥ⁺, Dᵥ⁻,T,mask)
    Λ = zeros(Float64,T^2,T^2)
    P = zeros(Float64,T^2)
    Q = zeros(Float64,T^2)
    Z⁰ = zeros(Float64,T^2)
    Z = zeros(Float64,T^2)
    # Vectorize inputs
    for j in CartesianIndices(z)
        Λ[j[1]+(j[2]-1)*T,j[1]+(j[2]-1)*T] = sqrt(λ[j])
        P[j[1]+(j[2]-1)*T] = p[j]
        Q[j[1]+(j[2]-1)*T] = q[j]
        Z⁰[j[1]+(j[2]-1)*T] = z⁰[j]
        Z[j[1]+(j[2]-1)*T] = z[j]
    end
    Λ = sparse(Diagonal(Λ))
    wᵤ⁺ = fill(1.0, T^2)
    wᵤ⁻ = fill(1.0, T^2)
    wᵥ⁺ = fill(1.0, T^2)
    wᵥ⁻ = fill(1.0, T^2)

    Wᵤ⁺ = Diagonal(wᵤ⁺)
    Wᵤ⁻ = Diagonal(wᵤ⁻)
    Wᵥ⁺ = Diagonal(wᵥ⁺)
    Wᵥ⁻ = Diagonal(wᵥ⁻)

    A1 = μ*(transpose(Dᵤ⁺)*transpose(Wᵤ⁺)*Wᵤ⁺*Dᵤ⁺ + transpose(Dᵤ⁻)*transpose(Wᵤ⁻)*Wᵤ⁻*Dᵤ⁻ + transpose(Dᵥ⁺)*transpose(Wᵥ⁺)*Wᵥ⁺*Dᵥ⁺ + transpose(Dᵥ⁻)*transpose(Wᵥ⁻)*Wᵥ⁻*Dᵥ⁻)
    A1 = A1 + Λ^2
    b = μ*((transpose(Dᵤ⁺)*transpose(Wᵤ⁺)*Wᵤ⁺)*P + (transpose(Dᵤ⁻)*transpose(Wᵤ⁻)*Wᵤ⁻)*P + (transpose(Dᵥ⁺)*transpose(Wᵥ⁺)*Wᵥ⁺)*Q + (transpose(Dᵥ⁻)*transpose(Wᵥ⁻)*Wᵥ⁻)*Q)
    b = b + (Λ^2)*Z⁰

    b₀ =  Diagonal(fill(1.0/(4.0*ϵ), T^2))
    DDᵤ⁺ = (ϵ*transpose(Dᵤ⁺)*Dᵤ⁺+b₀)
    DDᵤ⁻ = (ϵ*transpose(Dᵤ⁻)*Dᵤ⁻+b₀)
    DDᵥ⁺ = (ϵ*transpose(Dᵥ⁺)*Dᵥ⁺+b₀)
    DDᵥ⁻ = (ϵ*transpose(Dᵥ⁻)*Dᵥ⁻+b₀)
    A2 = ((Diagonal(Dᵤ⁻*Z-P)^2)*μ) + DDᵤ⁻
    b2 = fill(1.0/(4.0*ϵ), T^2)
    for i = 1:max_iter
        cg!(Z, A1, b)

        A2 = ((Diagonal(Dᵤ⁺*Z-P)^2)*μ) + DDᵤ⁺
        cg!(wᵤ⁺, A2, b2)
        Wᵤ⁺ = Diagonal(wᵤ⁺)

        A2 = ((Diagonal(Dᵤ⁻*Z-P)^2)*μ) + DDᵤ⁻
        cg!(wᵤ⁻, A2, b2)
        Wᵤ⁻ = Diagonal(wᵤ⁻)

        A2 = ((Diagonal(Dᵥ⁺*Z-Q)^2)*μ) + DDᵥ⁺
        cg!(wᵥ⁺, A2, b2)
        Wᵥ⁺ = Diagonal(wᵥ⁺)

        A2 = ((Diagonal(Dᵥ⁻*Z-Q)^2)*μ) + DDᵥ⁻
        cg!(wᵥ⁻, A2, b2)
        Wᵥ⁻ = Diagonal(wᵥ⁻)

        A1 = μ*(transpose(Dᵤ⁺)*transpose(Wᵤ⁺)*Wᵤ⁺*Dᵤ⁺ + transpose(Dᵤ⁻)*transpose(Wᵤ⁻)*Wᵤ⁻*Dᵤ⁻ + transpose(Dᵥ⁺)*transpose(Wᵥ⁺)*Wᵥ⁺*Dᵥ⁺ + transpose(Dᵥ⁻)*transpose(Wᵥ⁻)*Wᵥ⁻*Dᵥ⁻)
        A1 = A1 + Λ^2
        b = μ*((transpose(Dᵤ⁺)*transpose(Wᵤ⁺)*Wᵤ⁺)*P + (transpose(Dᵤ⁻)*transpose(Wᵤ⁻)*Wᵤ⁻)*P + (transpose(Dᵥ⁺)*transpose(Wᵥ⁺)*Wᵥ⁺)*Q + (transpose(Dᵥ⁻)*transpose(Wᵥ⁻)*Wᵥ⁻)*Q)
        b = b + (Λ^2)*Z⁰
    end

    for j in eachindex(index)
        z[index[j]] = Z[j]
    end
    z = z.*mask
    z = rotl90(z)
    return z
end



# function (scheme::Horn)(p::AbstractArray, q::AbstractArray, μ, ν)
#     @show "new"
#     iter = copy(scheme.max_iter)
#     ϵ = copy(scheme.ϵ)
#     Z = zeros(Float64, size(p))
#     Zᵏ⁺¹ = zeros(Float64, size(p))
#     R, C = size(p)
#     wᵤ⁺ = zeros(Float64, size(p))
#     wᵤ⁻ = zeros(Float64, size(p))
#     wᵥ⁺ = zeros(Float64, size(p))
#     wᵥ⁻ = zeros(Float64, size(p))
#     for i in CartesianIndices(p)
#         u,v=i.I
#         wᵤ⁺[i] = 1.0 / (sqrt((((q[min(u+1,R),v]-q[u,v])-(p[u,min(v+1,C)]-p[u,max(v-1,1)]))/μ)^2 + 1)*sqrt(((Z[min(u+1,R),v]-Z[u,v])^2/(ν^2)+1)))
#         wᵤ⁻[i] = 1.0 / (sqrt((((q[u,v]-q[max(u-1,1),v])-(p[u,min(v+1,C)]-p[u,max(v-1,1)]))/μ)^2 + 1)*sqrt(((Z[u,v]-Z[max(u-1,1),v])^2/(ν^2)+1)))
#         wᵥ⁺[i] = 1.0 / (sqrt((((q[min(u+1,R),v]-q[max(u-1,1),v])-(p[u,min(v+1,C)]-p[u,v]))/μ)^2 + 1)*sqrt(((Z[u,min(v+1,C)]-Z[u,v])^2/(ν^2)+1)))
#         wᵥ⁻[i] = 1.0 / (sqrt((((q[min(u+1,R),v]-q[max(u-1,1),v])-(p[u,v]-p[u,max(v-1,1)]))/μ)^2 + 1)*sqrt(((Z[u,v]-Z[u,max(v-1,1)])^2/(ν^2)+1)))
#     end
#     H = zeros(Float64, size(p))
#     V = zeros(Float64, size(p))
#     hv = zeros(Float64, size(p))
#     for i = 1:R
#         Z[i,1] = p[i,1]
#         Z[i,R] = p[i,R]
#     end
#     for i = 1:C
#         Z[1,i] = q[1,i]
#         Z[C,i] = q[C,i]
#     end
#     for i = 2:(R-1)
#         for j = 2:(C-1)
#             H[i,j] = (wᵥ⁻[i,j+1]*(p[i,j+1]-p[i,j]) + wᵥ⁺[i,j-1]*(p[i,j]-p[i,j-1])) / 2
#             V[i,j] = (wᵤ⁻[i+1,j]*(q[i+1,j]-q[i,j]) + wᵤ⁺[i-1,j]*(q[i,j]-q[i-1,j])) / 2
#             # h[i,j] = (p[i,j+1] - p[i,j-1]) / 2
#             # v[i,j] = (q[i+1,j] - q[i-1,j]) / 2
#             hv[i,j] = ϵ * (H[i,j] - V[i,j])
#             # hv[i,j] = hv[i,j]
#         end
#     end
#     for k = 1:iter
#         for i in CartesianIndices(p)
#             u,v=i.I
#             wᵤ⁺[i] = 1.0 / (sqrt((((q[min(u+1,R),v]-q[u,v])-(p[u,min(v+1,C)]-p[u,max(v-1,1)]))/μ)^2 + 1)*sqrt(((Z[min(u+1,R),v]-Z[u,v])^2/(ν^2)+1)))
#             wᵤ⁻[i] = 1.0 / (sqrt((((q[u,v]-q[max(u-1,1),v])-(p[u,min(v+1,C)]-p[u,max(v-1,1)]))/μ)^2 + 1)*sqrt(((Z[u,v]-Z[max(u-1,1),v])^2/(ν^2)+1)))
#             wᵥ⁺[i] = 1.0 / (sqrt((((q[min(u+1,R),v]-q[max(u-1,1),v])-(p[u,min(v+1,C)]-p[u,v]))/μ)^2 + 1)*sqrt(((Z[u,min(v+1,C)]-Z[u,v])^2/(ν^2)+1)))
#             wᵥ⁻[i] = 1.0 / (sqrt((((q[min(u+1,R),v]-q[max(u-1,1),v])-(p[u,v]-p[u,max(v-1,1)]))/μ)^2 + 1)*sqrt(((Z[u,v]-Z[u,max(v-1,1)])^2/(ν^2)+1)))
#         end
#         copyto!(Zᵏ⁺¹, Z)
#         for i = 2:(R-1)
#             for j = 2:(C-1)
#                 H[i,j] = (wᵥ⁻[i,j+1]*(p[i,j+1]-p[i,j]) + wᵥ⁺[i,j-1]*(p[i,j]-p[i,j-1])) / 2
#                 V[i,j] = (wᵤ⁻[i+1,j]*(q[i+1,j]-q[i,j]) + wᵤ⁺[i-1,j]*(q[i,j]-q[i-1,j])) / 2
#                 hv[i,j] = ϵ * (H[i,j] - V[i,j])
#                 Zᵏ⁺¹[i,j] = (wᵤ⁺[i-1,j]*Z[i-1,j] + wᵤ⁻[i+1,j]*Z[i+1,j] + wᵥ⁺[i,j-1]*Z[i,j-1] + wᵥ⁻[i,j+1]*Z[i,j+1])
#                 Zᵏ⁺¹[i,j] = Zᵏ⁺¹[i,j] - hv[i,j]
#                 Zᵏ⁺¹[i,j] = Zᵏ⁺¹[i,j] / (wᵤ⁻[i+1,j] + wᵤ⁺[i-1,j] + wᵥ⁻[i,j+1] + wᵥ⁺[i,j-1])
#             end
#         end
#         Z = Zᵏ⁺¹
#     end
#     return Z
# end
