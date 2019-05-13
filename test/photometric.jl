@testset "Photometric" begin
    I₁ = [0,0,1.0]
    I₂ = [1.0,0.2,1.0]
    I₃ = [0.2,0.9,1.0]
    img1, img2, img3 = generate_photometric(I₁, I₂, I₃, 1, shape=SynthSphere(), r = 5)
    Z, p, q = retrieve_surface(Photometric(), img1, img2, img3, -I₁, -I₂, -I₃)
    @test p[76,76] ≈ 0.25222893791413625
    @test p[1,1] ≈ 0.0
    @test p[151,151] ≈ 0.0
    @test q[90,90] ≈ 0.4603504763877097
    @test q[101,36] ≈ 0.930232558139535
    @test q[22, 136] ≈ 0.0
    @test Z[42,42] ≈ 3.879175028690482
end
