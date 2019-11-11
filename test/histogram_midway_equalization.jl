@testset "Midway Histogram Equalization" begin
    for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
        # Define two images consisting of luminance gradients with different
        # (non-overlapping) intensities. The "midway" equalization should produce
        # two new images for which the cummulative distribution function of the
        # image intensities is identical.
        img1 = T.(zeros(256,10))
        for i = 32:96
            for j = 1:10
                img1[i,j] = i/255.0
            end
        end
        img2 = T.(zeros(256,10))
        for i = 128:192
            for j = 1:10
                img2[i,j] = i/255.0
            end
        end
        img1o, img2o = adjust_histogram([img1, img2], MidwayEqualization(nbins = 256))
        edges1, counts1 = build_histogram(img1o, 256, minval = 0, maxval = 1)
        edges2, counts2 = build_histogram(img2o, 256, minval = 0, maxval = 1)
        @test sum(cumsum(counts2) - cumsum(counts1)) == 0

        edges, _ = build_histogram(img1, 256, minval = 0, maxval = 1)
        img1o, img2o = adjust_histogram([img1, img2], MidwayEqualization(edges = edges))
        edges1, counts1 = build_histogram(img1o, 256, minval = 0, maxval = 1)
        edges2, counts2 = build_histogram(img2o, 256, minval = 0, maxval = 1)
        @test sum(cumsum(counts2) - cumsum(counts1)) == 0

    end

    for T in (RGB{N0f8}, RGB{N0f16}, RGB{Float32}, RGB{Float64})
        imgg1 = zeros(Gray{Float32},256,10)
        for i = 32:96
            for j = 1:10
                imgg1[i,j] = i/255.0
            end
        end
        img1 = colorview(RGB,imgg1,imgg1,imgg1)
        img1 = T.(img1)

        imgg2 = zeros(Gray{Float32},256,10)
        for i = 128:192
            for j = 1:10
                imgg2[i,j] = i/255.0
            end
        end
        img2 = colorview(RGB,imgg2,imgg2,imgg2)
        img2 = T.(img2)

        img1o, img2o = ImageContrastAdjustment.adjust_histogram([img1, img2], ImageContrastAdjustment.MidwayEqualization(nbins = 256))

        edges1, counts1 = ImageContrastAdjustment.build_histogram(img1o, 256, minval = 0, maxval = 1)
        edges2, counts2 = ImageContrastAdjustment.build_histogram(img2o, 256, minval = 0, maxval = 1)
        # The algorithm equalizes the Y channels from the YIQ decomposition of the RGB images and then
        # constructs new RGB images by combining the equalised Y channels with the IQ components.
        # The build_histogram function then implicitly converts the "midway" RGB images to Gray
        # in order to construct the histogram. After this process the cummulative distribution functions
        # of these luminance gradients are no longer identical but still close.
        @test abs(sum(cumsum(counts2) - cumsum(counts1))) <= 20

        edges, _ = build_histogram(img1, 256, minval = 0, maxval = 1)
        img1o, img2o = adjust_histogram([img1, img2], MidwayEqualization(edges = edges))
        edges1, counts1 = build_histogram(img1o, 256, minval = 0, maxval = 1)
        edges2, counts2 = build_histogram(img2o, 256, minval = 0, maxval = 1)
        @test sum(cumsum(counts2) - cumsum(counts1)) == 0
    end
end
