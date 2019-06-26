@testset "pentland" begin
    img = generate_surface(SynthSphere(), 0.5, [0,0,1], radius = 50, background = false)
    Z = retrieve_surface(Pentland(), img)
    @test Z[76,76] ≈ 3.4232738028272176
    @test isapprox(Z[1,1], 1.160438047437746e-14; atol = 2e-15)
    @test isapprox(Z[151,151], 6.374017915783071e-15; atol = 2e-15)
    @test Z[90,90] ≈ 3.4196712241892793
    @test Z[101,36] ≈ 2.342218759464893
    @test isapprox(Z[22, 136], 9.89746106229774e-15; atol = 2e-15)
    @test Z[42,42] ≈ 2.9211132774577244
end
