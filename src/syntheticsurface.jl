function generate_surface(albedo::Real, illumination_direction::Vector{T} where T <: Real, radius::Real)
    ρ = albedo
    I = illumination_direction
    r = radius
    x,y = -1.5*r:0.1:1.5*r,-1.5*r:0.1:1.5*r
    if length(I) != 3
        @show "ADD ERROR CHECK HERE"
    end
    #p = -x./sqrt(r^2-(x.^2 + y.^2))
    #q = -y./sqrt(r^2-(x.^2 + y.^2))
    R = zeros(Complex{Float64}, length(x), length(y))
    for i = 1:length(x)
        for j = 1:length(y)
            p = -x[i]/sqrt(complex(r^2-(x[i]^2+y[j]^2)))
            q = -y[i]/sqrt(complex(r^2-(x[i]^2+y[j]^2)))
            R = (ρ*(-I[1]*p-I[2]*q+I[3]))/sqrt(1+p^2+q^2)
        end
    end
    @show sqrt(length(R))
    @show length(x)
end
