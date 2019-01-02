    @testset "Histogram Equalisation" begin

    for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
        #=
        Create an image that spans a narrow graylevel range. Then quantize
        the 256 bins down to 32 and determine how many bins have non-zero
        counts.
        =#

        img = Gray{Float32}.([i/255.0 for i = 64:128, j = 1:10])
        img = T.(img)
        _, counts_before = build_histogram(img, 32, minval = 0, maxval = 1)
        nonzero_before = sum(counts_before .!= 0)

        #=
        Equalize the image histogram. Then quantize the 256 bins down to 32
        and verify that all 32 bins have non-zero counts. This will confirm
        that the dynamic range of the original image has been increased.
        =#
        imgeq = adjust_histogram(Equalization(),img, 256, minval = 0, maxval = 1)
        edges, counts_after = build_histogram(imgeq, 32, minval = 0, maxval = 1)
        nonzero_after = sum(counts_after .!= 0)
        @test nonzero_before < nonzero_after
        @test nonzero_after == 32
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
        _, counts_before = build_histogram(img, 32, minval = 0, maxval = 1)
        nonzero_before = sum(counts_before .!= 0)

        #=
        Equalize the histogram. Then quantize the 256 bins down to 32 and
        verify that all 32 bins have non-zero counts. This will confirm that
        the dynamic range of the original image has been increased.
        =#
        imgeq = adjust_histogram(Equalization(),img,256, minval = 0, maxval = 1)
        edges, counts_after = build_histogram(imgeq, 32, minval = 0, maxval = 1)
        nonzero_after = sum(counts_after .!= 0)
        @test nonzero_before < nonzero_after
        @test nonzero_after == 32
    end

    # Verify that the minimum and maximum values of the equalised image match the
    # specified minimum and maximum values, i.e. that the intensities of the equalised
    # image are in the interval [minvalue, maxvalue].
    imgeq = adjust_histogram(Equalization(),collect(0:1:255), 256, minval = 64, maxval = 128)
    @test all(imgeq[1:65] .== 64)
    @test all(imgeq[128+1:end] .== 128)

    imgeq = adjust_histogram(Equalization(),collect(0:1/255:1), 256, minval = 64/255, maxval = 128/255)
    @test all(imgeq[1:65] .== 64/255)
    @test all(imgeq[128+1:end] .== 128/255)
end
