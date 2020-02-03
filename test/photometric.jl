@testset "Photometric" begin
    I₁ = [0,0,1.0]
    I₂ = [1.0,0.2,1.0]
    I₃ = [0.2,0.9,1.0]
    img1, img2, img3 = generate_photometric(SynthSphere(50), 1, I₁, I₂, I₃)
    Z, p, q = Photometric(-I₁, -I₂, -I₃, Horn())(img1, img2, img3)
    @test p[76,76] ≈ 0.25222893791413625
    @test p[1,1] ≈ 0.25222893791413625
    @test p[151,151] ≈ 0.25222893791413625
    @test q[90,90] ≈ 0.4603504763877097
    @test q[101,36] ≈ 0.930232558139535
    @test q[22, 136] ≈ 0.23815510035929432
    @test Z[42,42] ≈ 4.274056595786455
end
