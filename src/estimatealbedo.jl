function gradient(array::Array{T,2} where T <: Real)
    #calculate spatial gradient of array
    δx = similar(array)
    δy = similar(array)
    for i = 1:size(array,1)
        for j = 1:size(array,2)
            if i == 1
                δx[i,j] = array[i+1,j] - array[i,j]
            elseif i == size(array,1)
                δx[i,j] = array[i,j] - array[i-1,j]
            else
                δx[i,j] = (array[i+1,j] - array[i-1,j])/2
            end
            if j == 1
                δy[i,j] = array[i,j+1] - array[i,j]
            elseif j == size(array,2)
                δy[i,j] = array[i,j] - array[i,j-1]
            else
                δy[i,j] = (array[i,j+1] - array[i,j-1])/2
            end
        end
    end
    return δx, δy
end

function estimate_img_properties(img::AbstractArray)
    E = Array{Float64}(img)
    #calculate spatial gradient of img
    Ex,Ey = gradient(E)
    #normalize gradient
    nEx = similar(E)
    nEy = similar(E)
    Exy,~ = gradient(Ey)
    for i = 1:size(E,1)
        for j = 1:size(E,2)
            nEx[i,j] = Ex[i,j]/(Exy[i,j]+1^(-50))
            nEy[i,j] = Ey[i,j]/(Exy[i,j]+1^(-50))
        end
    end

    #calculate means
    μ₁ = mean(E)
    μ₂ = mean(E.^2)

    #calculate desired values
    g = sqrt(6*(π^2)*μ₂-48*(μ₁^2))
    ρ = g/π
    σ = acos((4*μ₁)/g)
    τ = atan(mean(nEy)/mean(nEx))
    if τ < 0
        τ = τ + π;
    end
    I = [cos(τ)*sin(σ),sin(τ)*sin(σ),cos(σ)]
    return ρ,I,σ,τ
end
