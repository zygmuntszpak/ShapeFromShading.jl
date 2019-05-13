function benchmark_iterative(Integration_Scheme::IntegrationScheme, iterations::AbstractArray, noiseLevels::AbstractArray, shape::SynthShape = SynthSphere(), r::Real = 5; logScale::Bool = true)
    Zₜ = ground_truth(shape, r)
    Zₜ = Zₜ ./ maximum(Zₜ).*35

    scene = Scene()
    runs = first(size(iterations))

    for i in eachindex(noiseLevels)
        p, q = sythetic_gradient(shape, radius = r, noise = noiseLevels[i])
        RMSE = zeros(Float64, 0)
        scene2 = Scene()
        Zout = zeros(Float64, size(p))
        for j in eachindex(iterations)
            println("Part $i: $j/$runs")
            Z = convert_gradient(Integration_Scheme, p, q, iterations[j])
            Z = Z./maximum(Z).*35
            Zout = Z
            E = 0.0
            for i in CartesianIndices(Z)
                E = E + (Zₜ[i] - Z[i])^2
            end
            push!(RMSE, sqrt(E / (first(size(Z)) * last(size(Z))))/(maximum(Zₜ) - minimum(Zₜ)))
        end

        iterationsDisplay = copy(iterations)
        if logScale == true
            iterationsDisplay = log10.(iterations)
        end

        range = 1.0:1:151
        scene2 = surface(range, range, Zout)

        c = RGB{Float64}(rand()%2,rand()%2,rand()%2)
        scatter!(iterationsDisplay, RMSE, markersize = 0.25, color = c)
        lines!(scene, iterationsDisplay, RMSE, color = c)
        display(vbox(scene, scene2))
        display(AbstractPlotting.PlotDisplay(), vbox(scene, scene2))
    end
end

function benchmark_noniterative(Integration_Scheme::IntegrationScheme, noiseLevels::AbstractArray, shape::SynthShape = SynthSphere(), r::Real = 5; trials::Int = 50)
    Zₜ = ground_truth(shape, r)
    Zₜ = Zₜ ./ maximum(Zₜ)

    RMSE = zeros(Float64, 0)
    errorBars = Array{Pair{Point{2,Float32},Point{2,Float32}}, 1}(undef, 0)
    scene = Scene()
    for i in eachindex(noiseLevels)
        Error = Array{Float64}(undef, 0)
        for j = 1:trials
            println("Part $i: $j/$trials")
            p, q = sythetic_gradient(shape, radius = r, noise = noiseLevels[i])
            Z = convert_gradient(Integration_Scheme,p,q)
            Z = Z./maximum(Z)
            E = 0.0
            for i in CartesianIndices(Z)
                E = E + (Zₜ[i] - Z[i])^2
            end
            push!(Error, sqrt(E / (first(size(Z)) * last(size(Z)))))
        end
        tempError = mean(Error)/(maximum(Zₜ) - minimum(Zₜ))
        push!(RMSE, tempError)
        push!(errorBars, Point2f0(noiseLevels[i], tempError - std(Error./(maximum(Zₜ) - minimum(Zₜ)))) => Point2f0(noiseLevels[i], tempError + std(Error./(maximum(Zₜ) - minimum(Zₜ)))))
    end
    scatter!(scene, noiseLevels, RMSE, markersize = 0.15, color = :black)
    linesegments!(scene, errorBars, color = :red)
    display(scene)
end

function compare_benchmark(Integration_Schemes::Array{IntegrationScheme,1}, noiseLevels::AbstractArray, shape = SynthSphere(), r::Real = 5; trails::Int = 50)
    Zₜ = ground_truth(shape, r)
    Zₜ = Zₜ./maximum(Zₜ)
    scene = Scene()
    runs = first(size(noiseLevels))

    for i in eachindex(Integration_Schemes)
        RMSE = zeros(Float64, 0)
        errorBars = Array{Pair{Point{2,Float32},Point{2,Float32}}, 1}(undef, 0)
        for j in eachindex(noiseLevels)
            Error = Array{Float64}(undef, 0)
            println("Part $i: $j/$runs")
            for k = 1:trails
                p, q = sythetic_gradient(shape, radius = r, noise = noiseLevels[j])
                Z = convert_gradient(Integration_Schemes[i],p,q)
                Z = Z./maximum(Z)
                E = 0.0
                for i in CartesianIndices(Z)
                    E = E + (Zₜ[i] - Z[i])^2
                end
                push!(Error, sqrt(E / (first(size(Z)) * last(size(Z)))))
            end
            tempError = mean(Error)/(maximum(Zₜ) - minimum(Zₜ))
            push!(RMSE, tempError)
            push!(errorBars, Point2f0(noiseLevels[j], tempError - std(Error./(maximum(Zₜ) - minimum(Zₜ)))) => Point2f0(noiseLevels[j], tempError + std(Error./(maximum(Zₜ) - minimum(Zₜ)))))
        end
        c = RGB{Float64}(rand()%2,rand()%2,rand()%2)
        scatter!(scene, noiseLevels, RMSE, markersize = 0.15, color = c)
        linesegments!(scene, errorBars, color = c)
        display(scene)
    end
end
