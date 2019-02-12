# Benchmark iterative integration scheme against ground truth.
function benchmark_iterative(Integration_Scheme::IntegrationScheme, iterations::AbstractArray, noiseLevels::AbstractArray, shape::SynthShape = SynthSphere(), r::Real = 50; logScale::Bool = true)
    # setup ground truth
    Zₜ = ground_truth(shape, radius = r)
    Zₜ = Zₜ ./ maximum(Zₜ).*35

    # setup plotting
    scene = Scene()
    runs = first(size(iterations))

    # For each noise level bechmark the scheme
    for i in eachindex(noiseLevels)

        # setup test
        p, q = sythetic_gradient(shape, radius = r, noise = noiseLevels[i])
        RMSE = zeros(Float64, 0)
        scene2 = Scene()
        Zout = zeros(Float64, size(p))

        for j in eachindex(iterations)
            println("Part $i: $j/$runs")

            # Run tests
            Z = convert_gradient(Integration_Scheme, p, q, iterations[j])
            Z = Z./maximum(Z).*35
            Zout = Z

            # calculate error
            E = 0.0
            for i in CartesianIndices(Z)
                E = E + (Zₜ[i] - Z[i])^2
            end
            push!(RMSE, sqrt(E / (first(size(Z)) * last(size(Z))))/(maximum(Zₜ) - minimum(Zₜ)))
        end

        # Change to log scale
        iterationsDisplay = copy(iterations)
        if logScale == true
            iterationsDisplay = log10.(iterations)
        end

        # Display reconstructed surface
        range = 1.0:1:151
        scene2 = surface(range, range, Zout)

        # Plot results
        c = RGB{Float64}(rand()%2, rand()%2, rand()%2)
        scatter!(iterationsDisplay, RMSE, markersize = 0.25, color = c)
        lines!(scene, iterationsDisplay, RMSE, color = c)
        display(Makie.vbox(scene, scene2))
        display(AbstractPlotting.PlotDisplay(), Makie.vbox(scene, scene2))
    end
end

# Benchmark non-iterative integration scheme against ground truth.
function benchmark_noniterative(Integration_Scheme::IntegrationScheme, noiseLevels::AbstractArray, shape::SynthShape = SynthSphere(), r::Real = 50; trials::Int = 50)
    # setup ground truth
    Zₜ = ground_truth(shape, radius = r)
    Zₜ = Zₜ ./ maximum(Zₜ)

    # setup plotting
    RMSE = zeros(Float64, 0)
    errorBars = Array{Pair{Point{2,Float32},Point{2,Float32}}, 1}(undef, 0)
    scene = Scene()

    # For each noise level bechmark the scheme
    for i in eachindex(noiseLevels)

        # Average across multiple trials
        Error = Array{Float64}(undef, 0)
        for j = 1:trials
            println("Part $i: $j/$trials")

            # Setup and Run tests
            p, q = sythetic_gradient(shape, radius = r, noise = noiseLevels[i])
            Z = convert_gradient(Integration_Scheme,p,q)
            Z = Z./maximum(Z)

            # Calculate error
            E = 0.0
            for i in CartesianIndices(Z)
                E = E + (Zₜ[i] - Z[i])^2
            end
            push!(Error, sqrt(E / (first(size(Z)) * last(size(Z)))))
        end

        # Average error and calculate error bars.
        tempError = mean(Error)/(maximum(Zₜ) - minimum(Zₜ))
        push!(RMSE, tempError)
        push!(errorBars, Point2f0(noiseLevels[i], tempError - std(Error./(maximum(Zₜ) - minimum(Zₜ)))) => Point2f0(noiseLevels[i], tempError + std(Error./(maximum(Zₜ) - minimum(Zₜ)))))
    end

    # Plot results
    scatter!(scene, noiseLevels, RMSE, markersize = 0.15, color = :black)
    linesegments!(scene, errorBars, color = :red)
    display(scene)
end

# Benchmark multiple integration schemes (using defult settings) against ground truth for comparison.
function compare_benchmark(Integration_Schemes::Array{IntegrationScheme,1}, noiseLevels::AbstractArray, shape = SynthSphere(), r::Real = 50; trails::Int = 50)
    # setup ground truth
    Zₜ = ground_truth(shape, radius = r)
    Zₜ = Zₜ./maximum(Zₜ)

    # setup plotting
    scene = Scene()
    runs = first(size(noiseLevels))

    # Benchmark each scheme
    for i in eachindex(Integration_Schemes)

        # Setup Error
        RMSE = zeros(Float64, 0)
        errorBars = Array{Pair{Point{2,Float32},Point{2,Float32}}, 1}(undef, 0)

        # For each noise level bechmark the scheme
        for j in eachindex(noiseLevels)
            Error = Array{Float64}(undef, 0)
            println("Part $i: $j/$runs")

            # Average across multiple trials
            for k = 1:trails
                # Setup and run test
                p, q = sythetic_gradient(shape, radius = r, noise = noiseLevels[j])
                Z = convert_gradient(Integration_Schemes[i],p,q)
                Z = Z./maximum(Z)

                # Calculate error
                E = 0.0
                for i in CartesianIndices(Z)
                    E = E + (Zₜ[i] - Z[i])^2
                end
                push!(Error, sqrt(E / (first(size(Z)) * last(size(Z)))))
            end

            # Average error and calculate error bars.
            tempError = mean(Error)/(maximum(Zₜ) - minimum(Zₜ))
            push!(RMSE, tempError)
            push!(errorBars, Point2f0(noiseLevels[j], tempError - std(Error./(maximum(Zₜ) - minimum(Zₜ)))) => Point2f0(noiseLevels[j], tempError + std(Error./(maximum(Zₜ) - minimum(Zₜ)))))
        end

        #plot result
        c = RGB{Float64}(rand()%2,rand()%2,rand()%2)
        scatter!(scene, noiseLevels, RMSE, markersize = 0.15, color = c)
        linesegments!(scene, errorBars, color = c)
        display(scene)
    end
end

# Benchmark reconsrtucted shape against ground truth and return error map
function benchmark_shape(Z::AbstractArray, ground::AbstractArray)
    # Setup test
    z = copy(Z)
    g = copy(ground)
    z = z./maximum(z)
    g = g./maximum(g)
    a = zeros(size(z))

    # Calculate error
    E = 0.0
    for i in CartesianIndices(Z)
        E = E + (g[i] - z[i])^2
        a[i] = sqrt((g[i] - z[i])^2)
    end
    error = sqrt(E / (first(size(Z)) * last(size(Z))))
    return error/(maximum(g) - minimum(g)), a./(maximum(g) - minimum(g))
end
