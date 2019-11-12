@testset "Contrast Stretching" begin

    for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
        #=
        Create an image that spans a narrow graylevel range.  Then
        quantize the 256 bins down to 16 and determine how many bins have
        non-zero counts.
        =#
        img = T.([i/255.0 for i = 64:128, j = 1:10])
        _, counts_before = build_histogram(img,16, minval = 0, maxval = 1)
        nonzero_before = sum(counts_before .!= 0)

        #=
        Stretch the histogram. Then quantize the 256 bins down to 16 and
        verify that all 16 bins have non-zero counts. This will confirm that
        the dynamic range of the original image has been increased.
        =#
        ret = adjust_histogram(img, ContrastStretching(t = 0.4, slope = 17))
        edges, counts_after = build_histogram(ret,16, minval = 0, maxval = 1)
        nonzero_after = sum(counts_after .!= 0)
        @test nonzero_before < nonzero_after
        @test nonzero_after == 16
        @test eltype(img) == eltype(ret)

        #=
        Verify that the function can cope with a NaN value.
        =#
        if T <: Gray{Float32} || T <: Gray{Float64}
            img[1] = NaN
            ret = adjust_histogram(img, ContrastStretching(t = 0.4, slope = 17))
            edges, counts_after = build_histogram(ret,16, minval = 0, maxval = 1)
            nonzero_after = sum(counts_after .!= 0)
            @test nonzero_before < nonzero_after
            @test nonzero_after == 16
            @test eltype(img) == eltype(ret)
        end

        #=
        Verify that when the slope is set to a very large value the contrast
        streching behaves like a thresholding function.
        =#
        ret = adjust_histogram(img, ContrastStretching(t = 0.37, slope = 1000))
        edges, counts_after = build_histogram(ret,16, minval = 0, maxval = 1)
        @test sum(counts_after .!= 0) == 2
    end

    for T in (RGB{N0f8}, RGB{N0f16}, RGB{Float32}, RGB{Float64})
        #=
        Create a color image that spans a narrow graylevel range.  Then
        quantize the 256 bins down to 16 and determine how many bins have
        non-zero counts.
        =#
        imgg = Gray{Float32}.([i/255.0 for i = 64:128, j = 1:10])
        img = colorview(RGB,imgg,imgg,imgg)
        img = T.(img)
        _, counts_before = build_histogram(img, 16, minval = 0, maxval = 1)
        nonzero_before = sum(counts_before .!= 0)

        #=
        Stretch the histogram. Then quantize the 256 bins down to 16 and
        verify that all 16 bins have non-zero counts. This will confirm that
        the dynamic range of the original image has been increased.
        =#
        ret = adjust_histogram(img, ContrastStretching(t = 0.4, slope = 17))
        edges, counts_after = build_histogram(ret, 16, minval = 0, maxval = 1)
        nonzero_after = sum(counts_after .!= 0)
        @test nonzero_before < nonzero_after
        @test nonzero_after == 16
        @test eltype(img) == eltype(ret)

        #=
        Verify that when the slope is set to a very large value the contrast
        streching behaves like a thresholding function.
        =#
        ret = adjust_histogram(img, ContrastStretching(t = 0.37, slope = 1000))
        edges, counts_after = build_histogram(ret, 16, minval = 0, maxval = 1)
        @test sum(counts_after .!= 0) == 2
    end
end
