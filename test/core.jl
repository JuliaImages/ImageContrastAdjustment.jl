@testset "Core" begin
        A = [NaN 1 2 3;
             NaN 6 5 4]
        @test ImageContrastAdjustment.minfinite(A) == 1
        @test ImageContrastAdjustment.maxfinite(A) == 6
        A = rand(10:20, 5, 5)
        @test ImageContrastAdjustment.minfinite(A) == minimum(A)
        @test ImageContrastAdjustment.maxfinite(A) == maximum(A)
        A = reinterpret(N0f8, rand(0x00:0xff, 5, 5))
        @test ImageContrastAdjustment.minfinite(A) == minimum(A)
        @test ImageContrastAdjustment.maxfinite(A) == maximum(A)
        A = rand(Float32,3,5,5)
        img = colorview(RGB, A)
        dc = ImageContrastAdjustment.minfinite(img)-RGB{Float32}(minimum(A, dims=(2,3))...)
        @test abs(dc) < 1e-6
        dc = ImageContrastAdjustment.maxfinite(img)-RGB{Float32}(maximum(A, dims=(2,3))...)
        @test abs(dc) < 1e-6
        @test ImageContrastAdjustment.minfinite(x->x^2,[NaN,10,2]) == 4
        @test ImageContrastAdjustment.maxfinite(x->x^2,[NaN,10,2]) == 100
end
