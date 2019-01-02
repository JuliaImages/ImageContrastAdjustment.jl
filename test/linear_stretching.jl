@testset "Linear Stretching" begin

    for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
        #=
        Stretching an image consisting of a linear ramp should not change the image
        if the specified minimum and maximum values match the minimum and maximum of
        the image.
        =#
        img = T.(collect(reshape(1/100:1/100:1, 10, 10)))
        minval = minimum(img)
        maxval = maximum(img)
        ret =  adjust_histogram(LinearStretching(), img, minval = minval, maxval = maxval)
        if T <: Gray{Float32} || T <: Gray{Float64}
            @test all(map((i, r) -> isapprox(i, r), img, ret))
        else
            @test ret == img
        end

        # Verify that NaN is also handled correctly.
        if T <: Gray{Float32} || T <: Gray{Float64}
            img[10] = NaN
            ret = adjust_histogram(LinearStretching(), img, minval = minval, maxval = maxval )
            @test isapprox(first(img), first(ret))
            @test isapprox(last(img), last(ret))
            @test isnan(ret[10])
        end

        # Verify that the smallest and largest values match the specified minval and maxval.
        img = T.(collect(reshape(1/100:1/100:1, 10, 10)))
        minval = minimum(img)
        maxval = maximum(img)
        ret =  adjust_histogram(LinearStretching(), img, minval = 0, maxval = 1)
        @test isapprox(0, first(ret))
        @test isapprox(1, last(ret))
        @test isapprox(0, minimum(ret[.!isnan.(ret)]))
        @test isapprox(1, maximum(ret[.!isnan.(ret)]))

        # Verify that the return type matches the input type.
        img = T.(testimage("mandril_gray"))
        ret = adjust_histogram(LinearStretching(), img, minval = 0, maxval = 1)
        @test eltype(ret) == eltype(img)
        @test isapprox(0, minimum(ret))
        @test isapprox(1, maximum(ret))

        ret = adjust_histogram(LinearStretching(), img, minval = 0.2, maxval = 0.8)
        @test eltype(ret) == eltype(img)
        @test isapprox(0.2, minimum(ret))
        @test isapprox(0.8, maximum(ret))
    end

    for T in (RGB{N0f8}, RGB{N0f16}, RGB{Float32}, RGB{Float64})
        #=
        Create a color image that spans a narrow graylevel range.  Then
        quantize the 256 bins down to 32 and determine how many bins have
        non-zero counts.
        =#
        imgg = Gray{Float32}.([i/255.0 for i = 64:128, j = 1:10])
        img = colorview(RGB,imgg,imgg,imgg)
        img = T.(img)
        _, counts_before = build_histogram(img,32, minval = 0, maxval = 1)
        nonzero_before = sum(counts_before .!= 0)

        #=
        Stretch the histogram. Then quantize the 256 bins down to 32 and
        verify that all 32 bins have non-zero counts. This will confirm that
        the dynamic range of the original image has been increased.
        =#
        ret = adjust_histogram(LinearStretching(),img, minval = 0, maxval = 1)
        edges, counts_after = build_histogram(ret, 32, minval = 0, maxval = 1)
        nonzero_after = sum(counts_after .!= 0)
        @test nonzero_before < nonzero_after
        @test nonzero_after == 32
        @test eltype(img) == eltype(ret)
    end

end
