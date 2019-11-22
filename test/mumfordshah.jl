@testset "MumfordShah" begin
    p, q = synthetic_gradient(SynthSphere(15), img_size = 50)
    z = zeros(Float64, size(p))
    Z = MumfordShah(z=z)(p, q)
    @test Z[25,25] ≈ 10.598521326817886
    @test Z[1,1] ≈ -2.2208982959969585
    @test Z[38,38] ≈ -2.1579445391253027
end
