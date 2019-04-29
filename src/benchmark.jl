function benchmark_iterative(Integration_Scheme::IntegrationScheme, iterations::AbstractArray, noiseLevels::AbstractArray, shape::SynthShape = SynthSphere(), r::Real = 5; logScale::Bool = true)
    Zₜ = ground_truth(shape, r)
    Zₜ = Zₜ ./ maximum(Zₜ) * 50

    scene = Scene()
    runs = first(size(iterations))

    for i in eachindex(noiseLevels)
        p, q = sythetic_gradient(shape, radius = r, noise = noiseLevels[i])
        RMSE = zeros(Float64, 0)

        for j in eachindex(iterations)
            println("Part $i: $j/$runs")
            Z = convert_gradient(Integration_Scheme, p, q, iterations[j])
            Z = Z./maximum(Z) * 50
            E = 0.0
            for i in CartesianIndices(Z)
                E = E + (Zₜ[i] - Z[i])^2
            end
            push!(RMSE, sqrt(E / (first(size(Z)) * last(size(Z)))))
        end

        iterationsDisplay = copy(iterations)
        if logScale == true
            iterationsDisplay = log10.(iterations)
        end

        c = RGB{Float64}(rand()%2,rand()%2,rand()%2)
        scatter!(iterationsDisplay, RMSE, markersize = 0.5, color = c)
        lines!(scene, iterationsDisplay, RMSE, color = c)
        display(scene)
    end
end

function benchmark_noniterative(Integration_Scheme::IntegrationScheme, noiseLevels::AbstractArray, shape::SynthShape = SynthSphere(), r::Real = 5; trails::Int = 50)
    Zₜ = ground_truth(shape, r)
    Zₜ = Zₜ ./ maximum(Zₜ) * 50

    RMSE = zeros(Float64, 0)
    scene = Scene()
    for i in eachindex(noiseLevels)
        Error = Array{Float64}(undef, 0)
        for j = 1:trails
            println("Part $i: $j/$trails")
            p, q = sythetic_gradient(shape, radius = r, noise = noiseLevels[i])
            Z = convert_gradient(Integration_Scheme,p,q)
            Z = Z./maximum(Z)*50
            E = 0.0
            for i in CartesianIndices(Z)
                E = E + (Zₜ[i] - Z[i])^2
            end
            push!(Error, sqrt(E / (first(size(Z)) * last(size(Z)))))
        end
        push!(RMSE, mean(Error))
    end
    scatter!(noiseLevels, RMSE, markersize = 0.5, color = :black)
    display(scene)
end

function compare_benchmark(Integration_Schemes::Array{IntegrationScheme,1}, noiseLevels::AbstractArray, shape = SynthSphere(), r::Real = 5; trails::Int = 50)
    Zₜ = ground_truth(shape, r)
    Zₜ = Zₜ./maximum(Zₜ)*50
    scene = Scene()
    runs = first(size(noiseLevels))

    for i in eachindex(Integration_Schemes)
        RMSE = zeros(Float64, 0)
        for j in eachindex(noiseLevels)
            Error = Array{Float64}(undef, 0)
            println("Part $i: $j/$runs")
            for k = 1:trails
                p, q = sythetic_gradient(shape, radius = r, noise = noiseLevels[j])
                Z = convert_gradient(Integration_Schemes[i],p,q)
                Z = Z./maximum(Z)*50
                E = 0.0
                for i in CartesianIndices(Z)
                    E = E + (Zₜ[i] - Z[i])^2
                end
                push!(Error, sqrt(E / (first(size(Z)) * last(size(Z)))))
            end
            push!(RMSE, mean(Error))
        end
        c = RGB{Float64}(rand()%2,rand()%2,rand()%2)
        scatter!(noiseLevels, RMSE, markersize = 0.5, color = c)
        display(scene)
    end
end
