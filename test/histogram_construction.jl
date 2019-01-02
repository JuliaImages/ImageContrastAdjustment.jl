@testset "Histogram Construction" begin
    # Consider an image where each intensity occurs only once and vary the number
    # of bins used in the histogram in powers of two. With the exception of the
    # first bin (with index 0), all other bins should have equal counts.
    expected_counts = [2^i for i = 0:7]
    bins = [2^i for i = 8:-1:1]
    for i = 1:length(bins)
        for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
            edges, counts  = build_histogram(T.(collect(0:1/255:1)),bins[i], minval = 0, maxval = 1)
            @test length(edges) == length(counts) - 1
            @test all(counts[1:end] .== expected_counts[i]) && counts[0] == 0
            @test axes(counts) == (0:length(edges),)
        end

        # Verify that the function can also take a color image as an input.
        for T in (RGB{N0f8}, RGB{N0f16}, RGB{Float32}, RGB{Float64})
            imgg = collect(0:1/255:1)
            img = colorview(RGB,imgg,imgg,imgg)
            edges, counts  = build_histogram(T.(img),bins[i], minval = 0, maxval = 1)
            @test length(edges) == length(counts) - 1
            @test all(counts[1:end] .== expected_counts[i]) && counts[0] == 0
            @test axes(counts) == (0:length(edges),)
        end

        # Consider also integer-valued images.
        edges, counts  = build_histogram(0:1:255,bins[i], minval = 0, maxval = 255)
        @test length(edges) == length(counts) - 1
        @test all(counts[1:end] .== expected_counts[i]) && counts[0] == 0
        @test axes(counts) == (0:length(edges),)
    end

    # Consider truncated intervals.
    for T in (Int, Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
        if T == Int
            edges, counts  = build_histogram(0:1:255,4, minval = 128, maxval = 192)
            @test length(edges) == length(counts) - 1
            @test collect(counts) == [128; 16; 16; 16; 80]
            @test axes(counts) == (0:length(edges),)
        else
            img = collect(0:1/255:1)
            edges, counts  = build_histogram(T.(img),4, minval = 128/255, maxval = 192/255)
            @test length(edges) == length(counts) - 1
            @test collect(counts) == [128; 16; 16; 16; 80]
            @test axes(counts) == (0:length(edges),)
        end

        if T == Int
            edges, counts  = build_histogram(0:1:255,4, minval = 120,maxval = 140)
            @test length(edges) == length(counts) - 1
            @test collect(counts) == [120, 5, 5, 5, 121]
            @test axes(counts) == (0:length(edges),)
        else
            img = collect(0:1/255:1)
            edges, counts  = build_histogram(T.(img),4,minval = 120/255, maxval = 140/255)
            @test length(edges) == length(counts) - 1
            @test axes(counts) == (0:length(edges),)
            # Due to roundoff errors the bins are not the same as in the
            # integer case above.
            if T == Gray{Float32}
                @test collect(counts) == [120, 5, 5, 5, 121]
            else
                @test collect(counts) == [120, 6, 5, 4, 121]
            end
        end
    end

    # Consider the case where the minimum and maximum values are not the start and
    # end points of the dynamic range. Because of numerical precision, the
    # results will be slightly different depending on the Image type.
    for T in (Int, Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
        if T == Int
            edges, counts  = build_histogram(200:1:240,4, minval = 200, maxval = 240)

            @test length(edges) == length(counts) - 1
            @test collect(counts) == [0, 10, 10, 10, 11]
            @test axes(counts) == (0:length(edges),)

            edges, counts  = build_histogram(200:1:240,4)
            @test length(edges) == length(counts) - 1
            @test collect(counts) == [0, 10, 10, 10, 11]
            @test axes(counts) == (0:length(edges),)
        else
            img = 200/255:1/255:240/255
            edges, counts  = build_histogram(T.(img),4, minval = 200/255, maxval = 240/255)
            @test length(edges) == length(counts) - 1
            if T == Gray{Float32}
                @test collect(counts) == [0, 10, 10, 10, 11]
            else
                @test collect(counts) == [0, 11, 9, 10, 11]
            end
            @test axes(counts) == (0:length(edges),)

            edges, counts  = build_histogram(T.(img),4)
            @test axes(counts) == (0:length(edges),)
            @test length(edges) == length(counts) - 1
            if T == Gray{N0f16} || T == Gray{N0f8}
                @test collect(counts) == [0, 11, 10, 10, 10]
            elseif T == Gray{Float32}
                @test collect(counts) == [0, 10, 10, 10, 11]
            else
                @test collect(counts) == [0, 11, 9, 10, 11]
            end

        end
    end

    # Consider the effect of NaN on the histogram counts.
    for j = 2:255
        for T in (Gray{Float32}, Gray{Float64})
            img = collect(0:1/255:1)
            img[j] = NaN
            edges, counts  = build_histogram(T.(img),256, minval = 0, maxval = 1)
            target = [1 for k = 1:length(counts)-1]
            target[j] -= 1
            @test length(edges) == length(counts) - 1
            @test counts[1:end] == target && counts[0] == 0
            @test axes(counts) == (0:length(edges),)

            # Verify that the minimum value (0) and maximum value (1) is
            # determined automatically even when there are NaN values
            # between these two bounds.
            img = collect(0:1/255:1)
            img[j] = NaN
            edges, counts  = build_histogram(T.(img),256)
            target = [1 for k = 1:length(counts)-1]
            target[j] -= 1 # Because of the NaN value in the image
            @test length(edges) == length(counts) - 1
            @test counts[1:end] == target && counts[0] == 0
            @test axes(counts) == (0:length(edges),)
        end
    end
end
