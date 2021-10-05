@testset "Piecewise Linear Stretching" begin

    @testset "constructors" begin
        img = testimage("mandril_gray")
        @test_throws ArgumentError PiecewiseLinearStretching((0..0.2 => 0..0.4, 0.2..0.2 => 0.4..1.0))
        @test_throws ArgumentError PiecewiseLinearStretching((0, 0.2, 0.2, 1), (0.0, 0, 0, 1))   
        @test_throws ArgumentError PiecewiseLinearStretching(src_knots = (0, 0.2, 1), dst_knots = (0.0, 0, 0, 1))
        @test_throws ArgumentError PiecewiseLinearStretching(src_knots = (0,), dst_knots = (1,))
        @test_throws ArgumentError PiecewiseLinearStretching(Percentiles(10,50, 100,150), (0, 0.5, 1), img)
    end

    @testset "miscellaneous" begin
        # We should be able to invert the image. 
         f = PiecewiseLinearStretching((0, 1), (1, 0))   
         for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
             img = T.(testimage("mandril_gray"))
             ret = adjust_histogram(img, f)
             @test ret ≈ (1 .- img)
         end 
        # We should be able to binarize the image. 
        f = PiecewiseLinearStretching((0.0..0.5 => 0..0, 0.5..1.0 => 1..1))
        for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
            img = T.(testimage("mandril_gray"))
            ret = adjust_histogram(img, f)
            vals = sort(unique(ret))
            @test first(vals) == 0.0
            @test last(vals) == 1.0
        end  

        # Setting the src_intervals equal to the dst_intervals should result in an identity
        # mapping (no change to the image) provided that the src_intervals span the entire
        # permissible contrast range, i.e. the unit interval. 
        f = PiecewiseLinearStretching((0, 0.3, 0.8 , 1.0), (0, 0.3, 0.8, 1.0))
        for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
            img = T.(testimage("mandril_gray"))
            ret = adjust_histogram(img, f)
            @test ret ≈  img
        end   

        # If the src_intervals equal the dst_intervals, but the src_intervals don't span the
        # entire permissible contrast range, i.e. the unit interval, then the resulting
        # image will not match the original image because pixels that fall outside the
        # specified src_intervals are saturated to the edge value of src_intervals if
        # saturate = true. 
        f = PiecewiseLinearStretching([0.3, 0.8], [0.3, 0.8]; saturate = true )
        for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
            img = T.(testimage("mandril_gray"))
            ret = adjust_histogram(img, f)
            @test !isapprox(ret, img, atol=1e-1)
        end  

        # If the src_intervals equal the dst_intervals, but the src_intervals don't span the
        # entire permissible contrast range, i.e. the unit interval, then the resulting
        # image will match the original image because pixels that fall outside the
        # specified src_intervals are not saturated to the edge value of src_intervals if
        # saturate = false
        f = PiecewiseLinearStretching([0.3, 0.8], [0.3, 0.8]; saturate = false )
        for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
            img = T.(testimage("mandril_gray"))
            ret = adjust_histogram(img, f)
            @test ret ≈ img
        end  

        for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})
            #=
            Stretching an image consisting of a linear ramp should not change the image
            if the specified minimum and maximum values match the minimum and maximum of
            the image.
            =#
            img = T.(collect(reshape(1/100:1/100:1, 10, 10)))
            minval = minimum(img)
            maxval = maximum(img)
            f = PiecewiseLinearStretching((minval..maxval => minval..maxval,))   
            ret =  adjust_histogram(img, f)
            if T <: Gray{Float32} || T <: Gray{Float64}
                @test all(map((i, r) -> isapprox(i, r), img, ret))
            else
                @test ret == img
            end

            # Verify that NaN is also handled correctly.
            if T <: Gray{Float32} || T <: Gray{Float64}
                img[10] = NaN
                ret = adjust_histogram(img, f)
                @test isapprox(first(img), first(ret))
                @test isapprox(last(img), last(ret))
                @test isnan(ret[10])
            end

            # Verify that the smallest and largest values match the specified minval and maxval.
            img = T.(collect(reshape(1/100:1/100:1, 10, 10)))
            minval = minimum(img)
            maxval = maximum(img)
            f = PiecewiseLinearStretching((minval..maxval => 0.0..1.0,)) 
            ret =  adjust_histogram(img, f)
            @test isapprox(0, first(ret))
            @test isapprox(1, last(ret))
            @test isapprox(0, minimum(ret[.!isnan.(ret)]))
            @test isapprox(1, maximum(ret[.!isnan.(ret)]))

            # Verify that the return type matches the input type.
            img = T.(testimage("mandril_gray"))
            minval = minimum(img)
            maxval = maximum(img)
            f = PiecewiseLinearStretching((minval..maxval => 0.0..1.0,)) 
            ret = adjust_histogram(img, f)
            @test eltype(ret) == eltype(img)
            @test isapprox(0, minimum(ret))
            @test isapprox(1, maximum(ret))

            f = PiecewiseLinearStretching((0.0..1.0 => 0.2..0.8,)) 
            ret = adjust_histogram(img, f)
            @test eltype(ret) == eltype(img)
            # We are mapping the interval [0, 1] to [0.2, 0.8] but since the 
            # smallest intensity in the mandril image need not be zero, and the largest
            # need to be one, smallest and largest observed values in the 
            # image need not be 0.2 and 0.8, but should still be within those bounds. 
            @test minimum(ret) >= 0.2 
            @test maximum(ret) <= 0.8 

            # Verify that results are correctly clamped to [0.2, 0.9] if it exceeds the range
            f = PiecewiseLinearStretching((0.3, 0.8), (0.2, 0.9))
            ret = adjust_histogram(img, f)
            @test eltype(ret) == eltype(img)
            @test minimum(ret) == T(0.2)
            @test maximum(ret) == T(0.9)  
        end

        # Verify that Percentiles and MinMax works correctly
        for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32}, Gray{Float64})   
            img = T.(collect(reshape(1/100:1/100:1, 10, 10)))
            f = PiecewiseLinearStretching(Percentiles((10,90)), MinMax(), img) 
            ret = adjust_histogram(img, f)
            @test all(isapprox.(ret[:,1], 0.01, atol = 1e-2))
            @test all(isapprox.(ret[:,end], 1.0, atol = 1e-2))

            # Same as above, but with keyword arguments.
            f = PiecewiseLinearStretching(img; src_knots = Percentiles((10,90)), dst_knots = MinMax()) 
            ret = adjust_histogram(img, f)
            @test all(isapprox.(ret[:,1], 0.01, atol=1e-2))
            @test all(isapprox.(ret[:,end], 1.0, atol=1e-2))

            f = PiecewiseLinearStretching(MinMax(), Percentiles((10,90)), img)
            ret = adjust_histogram(img, f)
            @test isapprox(ret[1,1], 0.109, atol=1e-2)
            @test isapprox(ret[end,end], 0.901, atol=1e-2)

            f = PiecewiseLinearStretching(Percentiles(10,90), Percentiles((10,90)), img)
            ret = adjust_histogram(img, f)
            @test all(isapprox.(ret[:,1], 0.109, atol=1e-2))
            @test all(isapprox.(ret[:,end], 0.901, atol=1e-2))
        end

        img = Float32.(collect(reshape(1/100:1/100:1, 10, 10)))

        f = PiecewiseLinearStretching((0.0,1.0), (-1.0, 0.0))
        ret = adjust_histogram(img, f)
        @test eltype(ret) == eltype(img)
        @test isapprox(minimum(ret), -1.0, atol=1e-1)
        @test isapprox(maximum(ret), 0.0, atol=1e-1)
        @test_throws DomainError adjust_histogram(Gray{N0f8}.(img), f)

        img = Gray{Float32}.(collect(reshape(1/100:1/100:1, 10, 10)))
        f = PiecewiseLinearStretching((0.0..0.2 => 0.0..0.2, 0.2..0.8 => -1.0..0.0, 0.8..1.0 => 1.0..(-2.0)))
        ret = adjust_histogram(img, f)
        @test eltype(ret) == eltype(img)
        @test minimum(ret) ≈ -2.0
        @test maximum(ret) ≈ 1.0

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
            f = PiecewiseLinearStretching((0.3, 0.5,), (0.0, 1.0))
            ret = adjust_histogram(img, f)
            edges, counts_after = build_histogram(ret, 32, minval = 0, maxval = 1)
            nonzero_after = sum(counts_after .!= 0)
            @test nonzero_before < nonzero_after
            @test nonzero_after == 32
            @test eltype(img) == eltype(ret)
        end    
    end
end