@testset "Gamma Correction" begin

    for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
        img = fill(oneunit(T),10,10)
        ret = adjust_histogram(GammaCorrection(), img, 1)
        @test img == ret
        @test eltype(ret) == eltype(img)

        imgp = padarray(img, Fill(zero(T), (2,2)))
        retp = adjust_histogram(GammaCorrection(), imgp, 1)
        @test imgp == retp
        @test eltype(retp) == eltype(imgp)
    end
    # ERROR: MethodError: no method matching ^(::AGray{Normed{UInt8,8}}, ::Float64)
    img = fill(oneunit(AGray{N0f8}),10,10)
    @test_broken adjust_histogram(GammaCorrection(),img, 0.5)

    # ERROR: MethodError: no method matching ^(::AGray{Normed{UInt8,8}}, ::Float64)
    img = fill(oneunit(ARGB{N0f8}),10,10)
    @test_broken adjust_histogram(GammaCorrection(),img, 0.5)


    for T in (RGB{N0f8}, RGB{N0f16}, RGB{Float64})
        img = fill(oneunit(T),10,10)
        ret = adjust_histogram(GammaCorrection(), img, 1)
        @test all(map((i, r) -> isapprox(i, r), img, ret))
        @test eltype(ret) == eltype(img)
        imgp = padarray(img, Fill(zero(T), (2,2)))
        retp = adjust_histogram(GammaCorrection(), imgp, 1)
        @test all(map((i, r) -> isapprox(i, r), imgp, retp))
        @test eltype(retp) == eltype(imgp)
    end

    for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
        img = T.(collect(reshape(1/100:1/100:1, 10, 10)))
        for i = 0.5:0.27:2
            ret = adjust_histogram(GammaCorrection(),img, i)
            @test ret == T.(img .^ i)
        end
    end

    # Since the function returns the same output type as the input type
    # there is an implicit rounding operation when dealing with integer values.
    img = reshape(1:1:100, 10, 10)
    for i = 0.5:0.28:2
        ret = adjust_histogram(GammaCorrection(),img, i)
        @test ret == round.(Int, img .^ i)
    end

end
