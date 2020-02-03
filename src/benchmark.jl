# Benchmark multiple integration schemes against ground truth for comparison.
function benchmark_integration(integration_schemes::Array{AbstractIntegrator,1}, noise_levels::AbstractArray, shape = SynthSphere(); trials::Int = 5, debug_mode = false)
    # setup ground truth
    Zₜ = ground_truth(shape)
    Zₜ = Zₜ./maximum(Zₜ)

    # setup plotting
    scene = Scene()
    runs = first(size(noise_levels))
    labels = string.(typeof.(integration_schemes))
    colours = Array{RGB{Float64}, 1}(undef, 0)
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
                p, q = synthetic_gradient(shape, noise = noise_levels[j])
                Z = integration_schemes[i](p,q)
                Z = Z./maximum(Z)

                # Calculate error
                E,a = benchmark_shape(Z, Zₜ)
                push!(error_list, E)
            end

            # Average error and calculate error bars.
            temp_error = mean(error_list)
            push!(RMSE, temp_error)
            push!(error_bars, Point2f0(noise_levels[j], temp_error - std(error_list./(maximum(Zₜ) - minimum(Zₜ)))) => Point2f0(noise_levels[j], temp_error + std(error_list./(maximum(Zₜ) - minimum(Zₜ)))))
        end

        #plot result
        c = RGB{Float64}(rand(),rand(),rand())
        lines!(scene, noise_levels, RMSE, color = c)
        # linesegments!(scene, error_bars, color = c)
        push!(colours , c)
        a = [scene[2]]
        if i > 2
            for k = 3:(i+1)
                push!(a, scene[k])
            end
        end
        lgd = legend(a, labels[1:i]; camera=campixel!, raw=true)
        s = Makie.vbox(scene, lgd)
        if !debug_mode
            display(s)
        end
    end
    a = [scene[2]]
    for i = 3:(length(labels)+1)
        push!(a, scene[i])
    end
    lgd = legend(a, labels; camera=campixel!, raw=true)
    scene = Makie.vbox(scene, lgd)
    if !debug_mode
        display(s)
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
    return error, a
end

function benchmark_image(img::AbstractArray, ground::AbstractArray)
    # Setup test
    R = Float64.(img)
    g = Float64.(ground)

    # Calculate error
    E = 0.0
    for i in CartesianIndices(R)
        E = E + (g[i] - R[i])^2
    end
    error = sqrt(E / (first(size(R)) * last(size(R))))
    return error
end

function benchmark_normals(P::AbstractArray, Q::AbstractArray, Pground::AbstractArray, Qground::AbstractArray)
    # Setup test
    p = copy(P)
    q = copy(Q)
    pground = copy(Pground)
    qground = copy(Qground)
    # Calculate error
    E = 0.0
    for i in CartesianIndices(p)
        E = E + (p[i] - pground[i])^2 + (q[i] - qground[i])^2
    end
    error = sqrt(E / (2 * first(size(p)) * last(size(q))))
    return error
end
