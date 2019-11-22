@testset "gen_mask" begin
    p, q = synthetic_gradient(Prism(), radius = 50)
    mask = gen_mask(p, q, 1)
    Z = Path()(p, q)
    @test mask[50, 50, 1] == 0.0
    @test mask[1, 1, 1] == 1.0
end
