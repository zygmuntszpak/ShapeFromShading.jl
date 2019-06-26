# Benchmark iterative integration scheme against ground truth.
function benchmark_iterative(integration_scheme::AbstractIntegrationScheme, iterations::AbstractArray, noise_levels::AbstractArray, shape::AbstractSyntheticShape = SynthSphere(), r::Real = 50; log_scale::Bool = true)
    # setup ground truth
    Zₜ = ground_truth(shape, radius = r)
    Zₜ = Zₜ ./ maximum(Zₜ).*35

    # setup plotting
    scene = Scene()
    runs = first(size(iterations))

    # For each noise level bechmark the scheme
    for i in eachindex(noise_levels)

        # setup test
        p, q = synthetic_gradient(shape, radius = r, noise = noise_levels[i])
        RMSE = zeros(Float64, 0)
        scene2 = Scene()
        Zout = zeros(Float64, size(p))

        for j in eachindex(iterations)
            println("Part $i: $j/$runs")

            # Run tests
            Z = convert_gradient(integration_scheme, p, q, iterations[j])
            Z = Z./maximum(Z).*35
            Zout = Z

            # Calculate error
            E = 0.0
            for i in CartesianIndices(Z)
                E = E + (Zₜ[i] - Z[i])^2
            end
            push!(RMSE, sqrt(E / (first(size(Z)) * last(size(Z))))/(maximum(Zₜ) - minimum(Zₜ)))
        end

        # Change to log scale
        iterations_display = copy(iterations)
        if log_scale == true
            iterations_display = log10.(iterations)
        end

        # Display reconstructed surface
        range = 1.0:1:151
        scene2 = surface(range, range, Zout)

        # Plot results
        c = RGB{Float64}(rand(), rand(), rand())
        scatter!(iterations_display, RMSE, markersize = 0.25, color = c)
        lines!(scene, iterations_display, RMSE, color = c)
        display(Makie.vbox(scene, scene2))
        display(AbstractPlotting.PlotDisplay(), Makie.vbox(scene, scene2))
    end
end

# Benchmark non-iterative integration scheme against ground truth.
function benchmark_noniterative(integration_scheme::AbstractIntegrationScheme, noise_levels::AbstractArray, shape::AbstractSyntheticShape = SynthSphere(), r::Real = 50; trials::Int = 50)
    # setup ground truth
    Zₜ = ground_truth(shape, radius = r)
    Zₜ = Zₜ ./ maximum(Zₜ)

    # setup plotting
    RMSE = zeros(Float64, 0)
    error_bars = Array{Pair{Point{2,Float32},Point{2,Float32}}, 1}(undef, 0)
    scene = Scene()

    # For each noise level bechmark the scheme
    for i in eachindex(noise_levels)

        # Average across multiple trials
        error_list = Array{Float64}(undef, 0)
        for j = 1:trials
            println("Part $i: $j/$trials")

            # Setup and Run tests
            p, q = synthetic_gradient(shape, radius = r, noise = noise_levels[i])
            Z = convert_gradient(integration_scheme,p,q)
            Z = Z./maximum(Z)

            # Calculate error
            E = 0.0
            for i in CartesianIndices(Z)
                E = E + (Zₜ[i] - Z[i])^2
            end
            push!(error_list, sqrt(E / (first(size(Z)) * last(size(Z)))))
        end

        # Average error and calculate error bars.
        temp_error = mean(error_list)/(maximum(Zₜ) - minimum(Zₜ))
        push!(RMSE, temp_error)
        push!(error_bars, Point2f0(noise_levels[i], temp_error - std(error_list./(maximum(Zₜ) - minimum(Zₜ)))) => Point2f0(noise_levels[i], temp_error + std(error_list./(maximum(Zₜ) - minimum(Zₜ)))))
    end

    # Plot results
    scatter!(scene, noise_levels, RMSE, markersize = 0.15, color = :black)
    linesegments!(scene, error_bars, color = :red)
    display(scene)
end

# Benchmark multiple integration schemes (using defult settings) against ground truth for comparison.
function compare_benchmark(integration_schemes::Array{AbstractIntegrationScheme,1}, noise_levels::AbstractArray, shape = SynthSphere(), r::Real = 50; trials::Int = 50)
    # setup ground truth
    Zₜ = ground_truth(shape, radius = r)
    Zₜ = Zₜ./maximum(Zₜ)

    # setup plotting
    scene = Scene()
    runs = first(size(noise_levels))

    # Benchmark each scheme
    for i in eachindex(integration_schemes)

        # Setup error_list
        RMSE = zeros(Float64, 0)
        error_bars = Array{Pair{Point{2,Float32},Point{2,Float32}}, 1}(undef, 0)

        # For each noise level bechmark the scheme
        for j in eachindex(noise_levels)
            error_list = Array{Float64}(undef, 0)
            println("Part $i: $j/$runs")

            # Average across multiple trials
            for k = 1:trials
                # Setup and run test
                p, q = synthetic_gradient(shape, radius = r, noise = noise_levels[j])
                Z = convert_gradient(integration_schemes[i],p,q)
                Z = Z./maximum(Z)

                # Calculate error
                E = 0.0
                for i in CartesianIndices(Z)
                    E = E + (Zₜ[i] - Z[i])^2
                end
                push!(error_list, sqrt(E / (first(size(Z)) * last(size(Z)))))
            end

            # Average error and calculate error bars.
            temp_error = mean(error_list)/(maximum(Zₜ) - minimum(Zₜ))
            push!(RMSE, temp_error)
            push!(error_bars, Point2f0(noise_levels[j], temp_error - std(error_list./(maximum(Zₜ) - minimum(Zₜ)))) => Point2f0(noise_levels[j], temp_error + std(error_list./(maximum(Zₜ) - minimum(Zₜ)))))
        end

        #plot result
        c = RGB{Float64}(rand(),rand(),rand())
        scatter!(scene, noise_levels, RMSE, markersize = 0.15, color = c)
        linesegments!(scene, error_bars, color = c)
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
