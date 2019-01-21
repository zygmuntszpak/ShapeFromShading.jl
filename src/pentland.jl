function retrieve_surface(algorithm::Pentland, img::AbstractArray)
    #find illumination and albedo
    ρ,I,σ,τ = estimate_img_properties(img)
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
