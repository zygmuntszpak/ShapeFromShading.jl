@testset "generate_surface" begin
    img = generate_surface(0.5, [0,0,1], radius = 5, scale_factor = 1.5, resolution=0.1)
    ρ,I,σ,τ = estimate_img_properties(img)
    @test ρ == 0.44074610935884634
    @test I ≈ [5.776768453507022e-17, 0.9434178830221143, 0.33160623937747646]
    @test σ == 1.2327906853082407
    @test τ == 1.5707963267948966
end