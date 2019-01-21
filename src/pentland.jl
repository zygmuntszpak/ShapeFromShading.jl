function retrieve_surface(algorithm::Pentland, img::AbstractArray)
    ρ,I,σ,τ = estimate_img_properties(img)
    return retrieve_surface(Pentland(),img,σ,τ)
end

#warning use at own risk
function retrieve_surface(algorithm::Pentland, img::AbstractArray, illumination_direction::Vector{T} where T <: Real)
    σ = acos(illumination_direction[3])
    τ = acos(illumination_direction[1]/sin(σ))
    @show σ,τ
    return retrieve_surface(Pentland(),img,σ,τ)
end

function retrieve_surface(algorithm::Pentland, img::AbstractArray, slant::Real, tilt::Real)
    #find illumination and albedo
    σ = slant
    τ = tilt
    E = Array{Float64}(img)
    #take fourier transform
    Fe = fft(E)
    M,N=size(E)

    #setup wx and wy
    wx,wy = setup_transform_values(M,N)

    #using the ilumination direction calculate the transformed Z
    Fz = Fe./(-1im.*wx.*cos(τ).*sin(σ)-1im.*wy.*sin(τ).*sin(σ))

    #recover Z
    Z = abs.(ifft(Fz))
    return Z
end
