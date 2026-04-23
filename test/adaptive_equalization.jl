@testset "Contrast Limited Adaptive Histogram Equalisation" begin
    for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
        # Create a graylevel ramp
        img = Gray{Float32}.([i/255.0 for i = 64:127, j = 1:64])
        img = T.(img)
        for rblocks in 2:4
            for cblocks in 2:4
                # There should be no change in intensity along a row.
                algo = AdaptiveEqualization(nbins = 256, minval = 0, maxval = 1,
                                            rblocks = rblocks, cblocks = cblocks,
                                            clip = 0.0)
                imgeq = adjust_histogram(img, algo)
                target = repeat(imgeq[:,1], inner=(1, 1), outer=(1, 64))
                @test all(imgeq .≈ target)

                # There should be no change in intensity along a column
                imgeq = adjust_histogram(transpose(img), algo)
                target =  transpose(repeat(imgeq[1, :], inner=(1, 1), outer=(1, 64)))
                @test all(imgeq .≈ target)
            end
        end

    end

    #=
    When rblocks and cblocks equal one and clip is zero,  then the method boils
    down to canonical histogram equalization. Hence, we repeat the same tests
    here that we have in histgoram_equalization.jl.
    =#

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
        algo = AdaptiveEqualization(nbins = 256, minval = 0, maxval = 1,
                                    rblocks = 1, cblocks = 1, clip = 0.0)
        imgeq = adjust_histogram(img, algo)
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
        algo = AdaptiveEqualization(nbins = 256, minval = 0, maxval = 1,
                                    rblocks = 1, cblocks = 1, clip = 0.0)
        imgeq = adjust_histogram(img, algo)
        edges, counts_after = build_histogram(imgeq, 32, minval = 0, maxval = 1)
        nonzero_after = sum(counts_after .!= 0)
        @test nonzero_before < nonzero_after
        @test nonzero_after == 32
    end

    # Verify that the minimum and maximum values of the equalised image match the
    # specified minimum and maximum values, i.e. that the intensities of the equalised
    # image are in the interval [minvalue, maxvalue].
    algo = AdaptiveEqualization(nbins = 256, minval = 64, maxval = 128,
                                rblocks = 1, cblocks = 1, clip = 0.0)
    imgeq = adjust_histogram(hcat(collect(0:1:255), collect(0:1:255)), algo)
    @test all(imgeq[1:65, :] .== 64)
    @test all(imgeq[128+1:end, :] .== 128)

    algo = AdaptiveEqualization(nbins = 256, minval = 64/255, maxval = 128/255,
                                rblocks = 1, cblocks = 1, clip = 0.0)
    imgeq = adjust_histogram(hcat(collect(0:1/255:1), collect(0:1/255:1)), algo)
    @test all(imgeq[1:65, :] .== 64/255)
    @test all(imgeq[128+1:end, :] .== 128/255)

    # Verify that increasing the clip value reduces the strength of the contrast
    # adjustment
    img = Gray{Float32}.(testimage("mandril_gray"))
    # As we increase the clip factor the contrast adjustment is diminished and we
    # should approach the original image.
    algo₁ = AdaptiveEqualization(nbins = 256, minval = 0, maxval = 1, rblocks = 1, cblocks = 1, clip = 0.0)
    algo₂ = AdaptiveEqualization(nbins = 256, minval = 0, maxval = 1, rblocks = 1, cblocks = 1, clip = 0.25)
    algo₃ = AdaptiveEqualization(nbins = 256, minval = 0, maxval = 1, rblocks = 1, cblocks = 1, clip = 0.5)
    algo₄ = AdaptiveEqualization(nbins = 256, minval = 0, maxval = 1, rblocks = 1, cblocks = 1, clip = 0.75)
    algo₅ = AdaptiveEqualization(nbins = 256, minval = 0, maxval = 1, rblocks = 1, cblocks = 1, clip = 1.0)
    algos = [algo₁, algo₂, algo₃, algo₄, algo₅]
    differences = [norm(adjust_histogram(img, algos[i]) .- img) for i = 1:5]
    @test all([differences[i] > differences[i+1] for i = 1:4])
    @test last(differences) < 2.8

    # Verify that auotmatically converting a color image and storing the histogram
    # adjusted result "inplace" works correctly.
    algo = AdaptiveEqualization(nbins = 256, minval = 0, maxval = 1, rblocks = 1, cblocks = 1, clip = 0)
    imgg₁ = Gray{Float32}.([i/255.0 for i = 64:128, j = 1:10])
    imgg₂ = copy(imgg₁)
    img = colorview(RGB,imgg₁,imgg₁,imgg₁)
    imgeq₁ = adjust_histogram!(imgg₁,img, algo)
    imgeq₂ = adjust_histogram(imgg₂, algo)
    @test norm(imgeq₁ .- imgeq₂) ≈ 0.0
end
