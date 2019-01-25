@testset "generate_surface" begin
    img = generate_surface(0.5, [0,0,1], radius = 5, scale_factor = 1.5, resolution=0.1)
    @test maximum(img) == 0.5
    @test size(img)[1] == 151
    @test img[10,10] == 0
end
