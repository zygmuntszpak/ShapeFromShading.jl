function generate_surface(albedo::Real = 0.5, illumination_direction::Vector{T} where T <: Real = [0, 0, 1]; radius::Real = 50, scale_factor::Real = 1.5, resolution::Real = 0.1)
    #initialize values
    ρ = albedo
    I = illumination_direction
    r = radius
    xyrange = -scale_factor*r:resolution:scale_factor*r
    x=zeros(length(xyrange),length(xyrange))
    y=zeros(length(xyrange),length(xyrange))
    for i in CartesianIndices(x)
        x[i]=xyrange[i[2]]
        y[i]=xyrange[i[1]]
    end
    R = zeros(Float64, axes(x))

    #calculate surface partial diferentials
    p = -x./sqrt.(complex(r^2 .-(x.^2+y.^2)))
    q = -y./sqrt.(complex(r^2 .-(x.^2+y.^2)))

    #calculate reflectance
    R = (ρ*(-I[1].*p-I[2].*q.+I[3]))./sqrt.(1 .+ p.^2+q.^2)

    #filter
    mask = ((r^2 .- (x.^2+y.^2) .> 0))
    R = R.*mask
    E = max.(0,Float64.(R))
    
    #convert to img and return
    img = colorview(Gray, E)
    return img
end
