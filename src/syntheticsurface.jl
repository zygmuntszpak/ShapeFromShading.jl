function generate_surface(albedo::Real = 0.5, illumination_direction::Vector{T} where T <: Real = [0, 0, 1]; radius::Real = 50, scale_factor::Real = 1.25, resolution::Real = 0.1)
    #initialize values
    ρ = albedo
    I = illumination_direction
    r = radius
    x = -scale_factor*r:resolution:scale_factor*r
    y= -scale_factor*r:resolution:scale_factor*r
    R = zeros(Float64, length(x), length(y))

    #create refection map
    for i = 1:length(x)
        for j = 1:length(y)
            #if point is inside circle continue else make = 0
            if r^2 - (x[i]^2 + y[j]^2) > 0
                #calculate surface partial derivatives at
                p = -x[i]/sqrt(r^2-(x[i]^2+y[j]^2))
                q = -y[j]/sqrt(r^2-(x[i]^2+y[j]^2))
                #calculate reflected value
                R[i,j] = (ρ*(-I[1]*p-I[2]*q+I[3]))/sqrt(1+p^2+q^2)
            else
                R[i,j] = 0.0
            end
            R[i,j] = max(0.0, R[i,j])
        end
    end

    #convert to img and return
    img = colorview(Gray, R)
    return img
end
