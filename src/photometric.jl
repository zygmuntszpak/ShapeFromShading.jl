function retrieve_surface(algorithm::Photometric, img1::AbstractArray, img2::AbstractArray, img3::AbstractArray, illumination_direction1::Vector{T} where T <: Real, illumination_direction2::Vector{T} where T <: Real, illumination_direction3::Vector{T} where T <: Real)
    #setup illumination matrix
    N = zeros(Float64, 3, 3)
    for i = 1:3
        N[1,i]=illumination_direction1[i]
        N[2,i]=illumination_direction2[i]
        N[3,i]=illumination_direction3[i]
    end
    Ninv = inv(N)

    #setup image and gradient arrays
    I₁ = Float64.(img1)
    I₂ = Float64.(img2)
    I₃ = Float64.(img3)
    p = zeros(Float64, size(I₁))
    q = zeros(Float64, size(I₁))

    #calculate gradients
    for i in CartesianIndices(I₁)
        I = [I₁[i], I₂[i], I₃[i]]
        n = Ninv * I
        if n[3] == 0.0
            n = [0,0,-1.0]
        end
        p[i] = -n[1] / n[3]
        q[i] = -n[2] / n[3]
    end

    #reconstruct surface
    Z = convert_gradient(p, q)
    # Z = mapwindow(median!, Z, (21,21))
    return Z, p, q
end
