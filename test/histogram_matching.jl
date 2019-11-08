@testset "Histogram Matching" begin
    # Construct a target image by adjusting the gamma of a standard
    # test image. Then transform the distribution of the test image
    # so that the intensity distribution matches the gamma-adjusted target
    # image. With Float32 and Float64 types the distributions match
    # exactly. For N0f8 and N0f16 types there are some round-off errors and
    # some of the bins of the transformed distributions are shifted
    # by one bin when compared to the target distribution. Hence, we
    # omit N0f8 and N0f16 from this particular test.
    for T in (Gray{Float32}, Gray{Float64})
        img = T.(testimage("mandril_gray"))
        imgo = T.(adjust_histogram(img, GammaCorrection(gamma = 0.5)))

        result = adjust_histogram(img, Matching(targetimg = imgo, nbins = 256))

        edges_target, counts_target = build_histogram(imgo, 256, minval = 0, maxval = 1)
        edges_result, counts_result = build_histogram(result, 256, minval = 0, maxval = 1)
        @test all(counts_target .== counts_result)

        # Check that we get the same result if we explicitly pass the
        # bin edges.
        edges, _ = build_histogram(img,256, minval = 0, maxval = 1)
        result = adjust_histogram(img, Matching(targetimg = imgo, edges = edges))

        edges_target, counts_target = build_histogram(imgo, 256,minval = 0, maxval = 1)
        edges_result, counts_result = build_histogram(result, 256, minval = 0, maxval = 1)
        @test all(counts_target .== counts_result)
    end

    for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float16}, Gray{Float32})
        img = T.([i < 64 ? i/255 : 0 for i = 0:255])
        imgo = T.([i > 64 && i < 128 ? i/255 : 0 for i = 0:255])
        result = adjust_histogram(img, Matching(targetimg = imgo, nbins = 256))
        edges_source, counts_source = build_histogram(img, 256, minval = 0, maxval = 1)
        edges_target, counts_target = build_histogram(imgo,256, minval = 0, maxval = 1)
        edges_result, counts_result = build_histogram(result,256, minval = 0, maxval = 1)
        @test all(counts_target .== counts_result)
    end

    for T in (RGB{N0f8}, RGB{N0f16}, RGB{Float32}, RGB{Float64})
        imgg = Gray{N0f8}.([i < 64 ? i/255 : 0 for i = 0:255])
        img =  T.(imgg,imgg,imgg)
        imggo = Gray{N0f8}.([i > 64 && i < 128 ? i/255 : 0 for i = 0:255])
        imgo = T.(imggo,imggo,imggo)
        result = adjust_histogram(img, Matching(targetimg = imgo, nbins = 256))
        edges_source, counts_source = build_histogram(img, 256, minval = 0, maxval = 1)
        edges_target, counts_target = build_histogram(imgo, 256, minval = 0, maxval = 1)
        edges_result, counts_result = build_histogram(result, 256, minval = 0, maxval = 1)
        @test all(counts_target .== counts_result)

        # Check that we get the same result if we explicitly pass the
        # bin edges.
        edges = 0.0:0.0019454658031463623:0.4960937798023224
        result = adjust_histogram(img, Matching(targetimg = imgo, edges = edges))
        edges_source, counts_source = build_histogram(img, 256, minval = 0, maxval = 1)
        edges_target, counts_target = build_histogram(imgo, 256, minval = 0, maxval = 1)
        edges_result, counts_result = build_histogram(result, 256, minval = 0, maxval = 1)
        @test all(counts_target .== counts_result)
    end

    # Modification of Tim's example:
    # https://github.com/JuliaImages/Images.jl/pull/752
    img = Gray.([0.4 0.45; 0.5 0.55])
    imgo = Gray.([0.1 0.3; 0.7 0.9])
    result = adjust_histogram(img, Matching(targetimg = imgo))
    edges_target, counts_target = build_histogram(imgo, 256, minval = 0, maxval = 1)
    edges_result, counts_result = build_histogram(result, 256, minval = 0, maxval = 1)
    @test all(counts_target .== counts_result)

    # Verify idempotency.
    img = collect(0:1:255)
    result = adjust_histogram(img, Matching(targetimg = img, nbins = 256))
    edges_target, counts_target = build_histogram(img, 256, minval = 0, maxval = 255)
    edges_result, counts_result = build_histogram(result, 256, minval = 0, maxval = 255)
    @test all(counts_target .== counts_result)

    # Verify that the algorithm works on integers.
    img = ([i <= 64 ? i : 0 for i = 0:255])
    imgo = ([i > 255-64 && i <= 255 ? i : 0 for i = 0:255])
    result = adjust_histogram(img, Matching(targetimg = imgo, nbins = 256))
    edges_source, counts_source = build_histogram(img, 256, minval = 0, maxval = 255)
    edges_target, counts_target = build_histogram(imgo, 256, minval = 0, maxval = 255)
    edges_result, counts_result = build_histogram(result, 256, minval = 0, maxval = 255)
    @test all(counts_target .== counts_result)

    # Verify that the algorithm can cope with NaN.
    img_nan =  Float64.(collect(0:1/255:1))
    img_nan[1:128] .= NaN
    reference_nan = adjust_histogram(img_nan, Matching(targetimg = img_nan))
    @test all(isapprox.(img_nan[129:end], reference_nan[129:end]; atol = 1e-1))
end
