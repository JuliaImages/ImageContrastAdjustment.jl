@testset "Core" begin
        A = [NaN 1 2 3;
             NaN 6 5 4]
        @test ImageContrastAdjustment.minimum_finite(A) == 1
        @test ImageContrastAdjustment.maximum_finite(A) == 6
        A = rand(10:20, 5, 5)
        @test ImageContrastAdjustment.minimum_finite(A) == minimum(A)
        @test ImageContrastAdjustment.maximum_finite(A) == maximum(A)
        A = reinterpret(N0f8, rand(0x00:0xff, 5, 5))
        @test ImageContrastAdjustment.minimum_finite(A) == minimum(A)
        @test ImageContrastAdjustment.maximum_finite(A) == maximum(A)
        A = rand(Float32,3,5,5)
        img = colorview(RGB, A)
        dc = ImageContrastAdjustment.minimum_finite(img)-RGB{Float32}(minimum(A, dims=(2,3))...)
        @test norm(dc) < 1e-6
        dc = ImageContrastAdjustment.maximum_finite(img)-RGB{Float32}(maximum(A, dims=(2,3))...)
        @test norm(dc) < 1e-6
        @test ImageContrastAdjustment.minimum_finite(x->x^2,[NaN,10,2]) == 4
        @test ImageContrastAdjustment.maximum_finite(x->x^2,[NaN,10,2]) == 100
end
