"""
```
    AdaptiveEqualization <: AbstractHistogramAdjustmentAlgorithm
    AdaptiveEqualization(; nbins = 256, minval = 0, maxval = 1, rblocks = 8, cblocks = 8, clip = 0.1)

    adjust_histogram([T,] img, f::AdaptiveEqualization)
    adjust_histogram!([out,] img, f::AdaptiveEqualization)
```

Performs Contrast Limited Adaptive Histogram Equalisation (CLAHE) on the input
image. It differs from ordinary histogram equalization in the respect that the
adaptive method computes several histograms, each corresponding to a distinct
section of the image, and uses them to redistribute the lightness values of the
image. It is therefore suitable for improving the local contrast and enhancing
the definitions of edges in each region of an image.

# Details

Histogram equalisation was initially conceived to  improve the contrast in a
single-channel grayscale image. The method transforms the distribution of the
intensities in an image so that they are as uniform as possible [1]. The natural
justification for uniformity is that the image has better contrast  if the
intensity levels of an image span a wide range on the intensity scale. As it
turns out, the necessary transformation is a mapping based on the cumulative
histogram---see [Equalization](@ref) for more details.

A natural extension of histogram equalisation is to apply the contrast
enhancement locally rather than globally [2]. Conceptually, one can imagine that
the process involves partitioning the image into a grid of rectangular regions
and applying histogram equalisation based on the local CDF of each contextual
region. However, to smooth the transition of the pixels from one contextual
region to another,  the mapping of a pixel is not necessarily done soley based
on the local CDF of its contextual region. Rather, the mapping of a pixel may be
interpolated based on the CDF of its contextual region, and the CDFs of the
immediate neighbouring regions.

In adaptive histogram equalisation the image ``\\mathbf{G}`` is partitioned into
``P \\times Q`` equisized submatrices,
```math
\\mathbf{G} =  \\begin{bmatrix}
\\mathbf{G}_{11} & \\mathbf{G}_{12} & \\ldots & \\mathbf{G}_{1C} \\\\
\\mathbf{G}_{21} & \\mathbf{G}_{22} & \\ldots & \\mathbf{G}_{2C} \\\\
\\vdots & \\vdots & \\ldots & \\vdots \\\\
\\mathbf{G}_{R1} & \\mathbf{G}_{R2} & \\ldots & \\mathbf{G}_{RC} \\\\
\\end{bmatrix}.
```

For each submatrix ``\\mathbf{G}_{rc}``, one computes a concomitant CDF, which we
shall denote by ``T_{rc}(G_{i,j})``. In the most general case, we will require
four CDFs
```math
\\begin{aligned}
T_1(v)  & \\triangleq  T_{rc}(G_{i,j}) \\\\
T_2(v)  & \\triangleq  T_{(r+1)c}(G_{i,j}) \\\\
T_3(v)  & \\triangleq  T_{(r+1)(c+1)}(G_{i,j}) \\\\
T_4(v)  & \\triangleq  T_{r(c+1)}(G_{i,j}).
\\end{aligned}
```

In order to determine which particular CDFs will be
used in the interpolation step, it is useful to (i) introduce the function

```math
\\Phi(\\mathbf{G}_{rc}) = \\left(  \\phi_{rc},  \\phi'_{rc}\\right) \\triangleq \\left(rP - \\frac{P}{2}, cQ - \\frac{Q}{2} \\right),
```

(ii) form the sequences  ``\\left(\\phi_{11}, \\phi_{21}, \\ldots, \\phi_{R1} \\right)``
and ``\\left(\\phi'_{11}, \\phi'_{12}, \\ldots, \\phi'_{1C} \\right)``, and (iii) define

```math
\\begin{aligned}
t  & \\triangleq  \\frac{i - \\phi_{r1}}{\\phi_{(r+1)1} - \\phi_{r1} } \\\\
u  & \\triangleq  \\frac{j - \\phi'_{1c}}{\\phi'_{1(c+1)} - \\phi'_{1c} }.
\\end{aligned}
```

#### Case I (Interior)
For a  pixel ``G_{i,j}`` in the range

```math
P - \\frac{P}{2} \\le i \\le RP - \\frac{P}{2}  \\quad \\text{and}  \\quad  Q - \\frac{Q}{2} \\le j \\le CQ - \\frac{Q}{2}.
```
values of ``r`` and ``c`` are implicitly defined by the solution to the inequalities
```math
\\phi_{r1} \\le i < \\phi_{(r+1)1}  \\quad \\text{and}  \\quad  \\phi'_{1c} \\le j < \\phi'_{1(c+1)}.
```

The *bilinearly interpolated* transformation that maps an intensity ``v`` at location ``(i,j)`` in the image
to an intensity ``v'`` is given by [3]

```math
v' \\triangleq \\bar{T}(v)  = (1-t) (1-u)T_1(G_{i,j}) + t(1-u)T_2(G_{i,j}) + tuT_3(G_{i,j}) +(1-t)uT_4(G_{i,j}).
```

#### Case II (Vertical Border)

For a  pixel ``G_{i,j}`` in the range

```math
P - \\frac{P}{2} \\le i \\le RP - \\frac{P}{2}  \\quad \\text{and}  \\quad  1 \\le j < Q - \\frac{Q}{2}  \\; \\cup \\;   CQ - \\frac{Q}{2} < j \\le CQ,
```

``r`` is implicitly defined by the solution to the inequality ``\\phi_{r1} \\le i < \\phi_{(r+1)1}``, while

```math
c = \\begin{cases}
   1, & \\text{if }  \\quad  1 \\le j < Q - \\frac{Q}{2}  \\\\
   C, & \\text{if } \\quad   CQ - \\frac{Q}{2} < j \\le CQ.
\\end{cases}
```

The *linearly interpolated* transformation that maps an intensity ``v`` at location ``(i,j)`` in the image
to an intensity ``v'`` is given by

```math
v' \\triangleq \\bar{T}(v)  = (1-t)T_1(G_{i,j}) + tT_2(G_{i,j}).
```

#### Case III (Horizontal Border)

For a  pixel ``G_{i,j}`` in the range

```math
1 \\le i < P - \\frac{P}{2}  \\;\\cup \\;   RP - \\frac{P}{2} < i \\le RP    \\quad \\text{and}  \\quad  Q - \\frac{Q}{2} \\le j \\le CQ - \\frac{Q}{2},
```

``c`` is implicitly defined by the solution to the inequality ``\\phi'_{1c} \\le j < \\phi'_{1(c+1)}``, while
```math
r = \\begin{cases}
   1, & \\text{if }  \\quad  1 \\le i < P - \\frac{P}{2}  \\\\
   R, & \\text{if } \\quad   RP - \\frac{P}{2} < i \\le RP .
\\end{cases}
```

The *linearly interpolated* transformation that maps an intensity ``v`` at location ``(i,j)`` in the image
to an intensity ``v'`` is given by

```math
v' \\triangleq \\bar{T}(v)  = (1-u)T_1(G_{i,j}) + uT_4(G_{i,j}).
```

#### Case IV (Corners)
For a  pixel ``G_{i,j}`` in the range

```math
1 \\le i < \\frac{P}{2}  \\;\\cup \\; RP - \\frac{P}{2} < i \\le RP   \\quad \\text{and}  \\quad  1 \\le j < CQ -  \\frac{Q}{2} \\; \\cup \\;   CQ - \\frac{Q}{2} < j \\le CQ ,
```
we have
```math
r = \\begin{cases}
   1, & \\text{if }  \\quad  1 \\le i < P - \\frac{P}{2}  \\\\
   R, & \\text{if } \\quad   RP - \\frac{P}{2} < i \\le RP
\\end{cases}
 \\quad \\text{and}  \\quad
c = \\begin{cases}
   1, & \\text{if }  \\quad  1 \\le j < Q - \\frac{Q}{2}  \\\\
   C, & \\text{if } \\quad   CQ - \\frac{Q}{2} < j \\le CQ.
\\end{cases}
```
The transformation that maps an intensity ``v`` at location ``(i,j)`` in the image
to an intensity ``v'`` is given by

```math
v' \\triangleq \\bar{T}(v)  = T_1(G_{i,j}).
```

## Limiting Contrast
An unfortunate side-effect of contrast enhancement is that it has a tendency to
amplify the level of noise in an image, especially when the magnitude of the
contrast enhancement is very high. The magnitude of contrast enhancement is
associated with the gradient of ``T(\\cdot)``, because the  gradient determines the
extent to which consecutive input intensities are stretched across the
grey-level spectrum. One can diminish the level of noise amplification by
limiting the magnitude of the contrast enhancement, that is, by limiting the
magnitude of the gradient.

Since the derivative of ``T(\\cdot)`` is the empirical density ``\\hat{f}_{G}``,
the slope of the mapping function at any input intensity is proportional to the
height of the histogram  ``\\hat{f}_{G}`` at that intensity.  Therefore,
limiting the slope of the local mapping function is equivalent to clipping the
height of the histogram. A detailed description of the  implementation  details
of the clipping process can be found in [2].

# Options

Various options for the parameters of this function are described in more detail
below.

## Choices for `img`

The function can handle a variety of input types. The returned image
depends on the input type.

For coloured images, the input is converted to
[YIQ](https://en.wikipedia.org/wiki/YIQ) type and the Y channel is equalised.
This is the combined with the I and Q channels and the resulting image converted
to the same type as the input.

## Choices for `nbins` in `AdaptiveEqualization`

You can specify the total number of bins in the histogram of each local region.

## Choices for `rblocks` and `cblocks` in `AdaptiveEqualization`

The `rblocks` and `cblocks` specify the number of blocks to divide the input
image into in each direction. By default both values are set to eight.

## Choices for `clip` in `AdaptiveEqualization`

The `clip` parameter must be a value between 0 and 1. It defines an implicit
threshold at which a histogram is clipped. Counts that exceed the threshold are
redistributed as equally as possible so that no bin exceeds the threshold limit.
A value of zero means no clipping, whereas a value of one sets the threshold at
the smallest feasible bin limit. A bin limit is feasible if all bin counts can be
redistributed such that no bin count exceeds the limit. In practice, a `clip` value
of zero corresponds to maximal contrast enhancement, whereas a `clip` value of
one corredponds to minimal contrast enhancement. The default value is `0.1`.

## Choices for `minval` and `maxval` in `AdaptiveEqualization`

If `minval` and `maxval` are specified then intensities are equalized to the range
[`minval`, `maxval`]. The default values are 0 and 1.

# Example

```julia

using TestImages, FileIO, ImageView

img =  testimage("mandril_gray")
imgeq = adjust_histogram(img, AdaptiveEqualization(nbins = 256, rblocks = 4, cblocks = 4, clip = 0.2))

imshow(img)
imshow(imgeq)
```

# References
1. R. C. Gonzalez and R. E. Woods. *Digital Image Processing (3rd Edition)*.  Upper Saddle River, NJ, USA: Prentice-Hall,  2006.
2. S. M. Pizer, E. P. Amburn, J. D. Austin, R. Cromartie, A. Geselowitz, T. Greer, B. ter Haar Romeny, J. B. Zimmerman and K. Zuiderveld “Adaptive histogram equalization and its variations,” *Computer Vision, Graphics, and Image Processing*, vol. 38, no. 1, p. 99, Apr. 1987. [10.1016/S0734-189X(87)80186-X](https://doi.org/10.1016/s0734-189x(87)80156-1)
3. W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery.  *Numerical Recipes: The Art of Scientific Computing (3rd Edition)*. New York, NY, USA: Cambridge University Press, 2007.
"""
@with_kw struct AdaptiveEqualization{T₁ <: Union{Real,AbstractGray},
                                        T₂ <: Union{Real,AbstractGray},
                                        T₃ <: Real} <: AbstractHistogramAdjustmentAlgorithm
    nbins::Int = 256
    minval::T₁ = 0.0
    maxval::T₂ = 1.0
    rblocks::Int = 8
    cblocks::Int = 8
    clip::T₃ = 0.1
end

function (f::AdaptiveEqualization)(out::GenericGrayImage, img::GenericGrayImage)
    validate_parameters(f)
    height, width = length.(axes(img))
    # If necessary, resize the image so that the requested number of blocks fit exactly.
    resized_height = ceil(Int, height / (2 * f.rblocks)) * 2 * f.rblocks
    resized_width = ceil(Int, width / (2 * f.cblocks)) * 2 * f.cblocks
    must_resize = (resized_height != height) || (resized_width != width) ? true : false
    if must_resize
        img_tmp = imresize(img, (resized_height, resized_width))
        out_tmp = copy(img_tmp)
    else
        img_tmp = img
        out_tmp = out
    end

    # Determine the relevant pixel coordinates for each block.
    block_height = resized_height ÷  f.rblocks
    block_width =  resized_width ÷  f.cblocks
    r_intervals = [StepRange((r - 1) * block_height + 1, 1, r * block_height) for r = 1:f.rblocks]
    c_intervals = [StepRange((c - 1) * block_width + 1, 1, c * block_width) for c = 1:f.cblocks]

    function construct_hist(roi::AbstractArray)
        edges, histogram = build_histogram(roi, f.nbins, minval = f.minval, maxval = f.maxval)
    end

    function construct_cdf(edges::AbstractArray, histogram::AbstractArray)
        lb = first(axes(histogram, 1))
        ub = last(axes(histogram, 1))
        csum = cumsum(histogram[lb:ub])
        cdf = csum / csum[end]
        return edges, cdf
    end

    # Construct a histogram for each block in the image.
    block_hist = [construct_hist(view(img_tmp, r_intervals[r], c_intervals[c], :)) for r = 1:f.rblocks, c = 1:f.cblocks]

    # Redistribute histogram counts in accordance with the clip weight.
    map!(x-> (x[1], clip_histogram!(x[2], f.clip)), block_hist, block_hist)

    # Construct a CDF for each block in the image.
    block_cdf = [construct_cdf(block_hist[r,c]...) for r = 1:f.rblocks, c = 1:f.cblocks]

    block_centroid_r = [r * block_height - block_height ÷ 2 for r = 1:f.rblocks]
    block_centroid_c = [c * block_width - block_width ÷ 2 for c = 1:f.cblocks]
    intensity_range = (f.minval, f.maxval)
    # Transform pixels using linear and bilinear interpolation of the block CDFs.
    transform_image!(out_tmp, img_tmp, block_centroid_r, block_centroid_c,
                     block_width, block_height, intensity_range,
                     f.cblocks, f.rblocks, block_cdf)
    out .= must_resize ?  imresize(out_tmp, (height, width)) : out_tmp
    return out
end

function validate_parameters(f::AdaptiveEqualization)
    !(0 <= f.clip <= 1) && throw(ArgumentError("The parameter `clip` must be in the range [0..1]."))
    !(1 <= f.rblocks && 1 <= f.cblocks) && throw(ArgumentError("The parameters `rblocks` and `cblocks` must be greater than 0."))
end

function (f::AdaptiveEqualization)(out::AbstractArray{<:Color3}, img::AbstractArray{<:Color3})
    T = eltype(img)
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    #=
       TODO: Understand the cause and solution of this error.
       When I pass a view I run into this error on Julia 1.1.
       ERROR: ArgumentError: an array of type `Base.ReinterpretArray` shares memory with another argument and must
       make a preventative copy of itself in order to maintain consistent semantics,
       but `copy(A)` returns a new array of type `Array{Float64,3}`. To fix, implement:
       `Base.unaliascopy(A::Base.ReinterpretArray)::typeof(A)`
    =#
    #adjust_histogram!(view(yiq_view,1,:,:), f)
    y = comp1.(yiq)
    adjust_histogram!(y, f)
    yiq_view[1, :, :] .= y
    out .= convert.(T, yiq)
end



(f::AdaptiveEqualization)(out::GenericGrayImage, img::AbstractArray{<:Color3}) =
    f(out, of_eltype(Gray, img))

function clip_histogram!(histogram::AbstractArray, clip_weight::Number)
    limit, unfilled_bin_count = determine_threshold(histogram, clip_weight)
    initial_excess = determine_excess!(histogram, limit)
    remaining_excess = perform_initial_redistribution!(histogram, limit, initial_excess, unfilled_bin_count)
    perform_iterative_redistribution!(histogram, limit, remaining_excess)
    return histogram
end

function determine_excess!(histogram::AbstractArray, limit::Number)
    excess = zero(limit)
    for n in eachindex(histogram)
        val = histogram[n]
        if val > limit
            excess = excess + (val - limit)
            histogram[n] = limit
        end
    end
    return excess
end

function determine_threshold(histogram::AbstractArray, clip_weight::Number)
    sorted_hist = sort(collect(histogram), rev = true)
    required_capacity = zeros(axes(sorted_hist))
    available_capacity = zeros(axes(sorted_hist))
    # If we pick an empirical bin count as the actual limit, determine how many
    # counts we will need to move (required_capacity) and how much space do
    # we have in the bins (available_capacity).
    for i = 1:length(sorted_hist)
        required_capacity[i] = sum(sorted_hist[1:i] .- sorted_hist[i])
        available_capacity[i] = sum(sorted_hist[i] .- sorted_hist[i:end])
    end
    # Find the smallest feasible bin limit.
    largest_feasible_index = 1
    for i in eachindex(required_capacity)
        if available_capacity[i] < required_capacity[i]
            largest_feasible_index = i - 1
            break
        end
    end
    # The target bin limit is a convex combination of the highest bin count
    # and the smallest feasible bin count.
    smallest_feasible_limit = sorted_hist[largest_feasible_index]
    target_limit = (1 - clip_weight) * first(sorted_hist) + (clip_weight) * smallest_feasible_limit
    # Set the limit to be the empirical count that is close to the target limit.
    chosen_index = searchsortedlast(sorted_hist, round(Int,target_limit); rev = true)
    threshold = sorted_hist[chosen_index]
    # Knowing how many bins have spare capacity will help us distribute the
    # excess bin counts more efficiently in the perform_initial_redistribution!
    # function.
    unfilled_bin_count  = length(sorted_hist) - chosen_index
    return threshold, unfilled_bin_count
end

function perform_initial_redistribution!(histogram::AbstractArray, limit::Number, excess::Number, N::Integer)
    m = excess ÷ N
    for n in eachindex(histogram)
        val = histogram[n]
        if excess > 0
            if val < limit - m
                histogram[n] = histogram[n] + m
                excess = excess - m
            elseif val < limit
                excess = excess - (limit - val)
                histogram[n] = limit
            end
        end
    end
    return excess
end

function perform_iterative_redistribution!(histogram::AbstractArray, limit::Number, excess::Number)
    while excess > 0
        for n in eachindex(histogram)
            val = histogram[n]
            if excess > 0
                if val < limit
                    excess = excess - 1
                    histogram[n] = histogram[n]  + 1
                end
            end
        end
    end
end


function apply_cdf_transform(val::Union{Real,AbstractGray}, minval::Union{Real,AbstractGray}, maxval::Union{Real,AbstractGray}, edges::AbstractArray, cdf::AbstractArray)
    first_edge = first(edges)
    inv_step_size = 1 / step(edges)
    scale = (maxval - minval) / (cdf[end] - first(cdf))
    if val >= edges[end]
        newval = cdf[end]
    elseif val < first_edge
        newval = first(cdf)
    else
        index = floor(Int, (val - first_edge) * inv_step_size) + 1
        newval = cdf[index]
    end
    # Scale the new intensity value to so that it lies in the range [minval, maxval].
    newval = minval + (newval - first(cdf)) * scale
end

function transform_image!(out, img, block_centroid_r, block_centroid_c, block_width, block_height, intensity_range, cblocks, rblocks, block_cdf)
    height, width = length.(axes(out))
    r₀ = first(block_centroid_r) + 1
    r₁ = last(block_centroid_r) - 1
    c₀ = first(block_centroid_c) + 1
    c₁ = last(block_centroid_c) - 1
    block_dimensions = (block_width, block_height)

    # Interior
    bounds = (r₀:r₁, c₀:c₁)
    block_centroids = (block_centroid_r, block_centroid_c)
    transform_interior!(out, img, bounds, block_centroids, block_dimensions,
                        intensity_range, block_cdf)
    # West
    bounds = (r₀:r₁, 1:(c₀-1))
    block_c_idx = 1
    transform_vertical_strip!(out, img, bounds, block_centroid_r, block_c_idx,
                              block_height, intensity_range, block_cdf)
    # East
    bounds = (r₀:r₁, (c₁+1):width)
    block_c_idx = cblocks
    transform_vertical_strip!(out, img, bounds, block_centroid_r, block_c_idx,
                              block_height, intensity_range, block_cdf)
    # North
    bounds = (1:r₀-1, c₀:c₁)
    block_r_idx = 1
    transform_horizontal_strip!(out, img, bounds, block_centroid_c, block_r_idx,
                                block_width, intensity_range, block_cdf)
    # South
    bounds = (r₁+1:height, c₀:c₁)
    block_r_idx = rblocks
    transform_horizontal_strip!(out, img, bounds, block_centroid_c, block_r_idx,
                                block_width, intensity_range, block_cdf)
    # North-West
    bounds = (1:(r₀ - 1), 1:(c₀ - 1))
    block_r_idx = 1
    block_c_idx = 1
    transform_corner!(out, img, bounds, block_r_idx, block_c_idx,
                      intensity_range, block_cdf)
    # North-East
    bounds = (1:(r₀ - 1), (c₁ + 1):width)
    block_r_idx = 1
    block_c_idx = cblocks
    transform_corner!(out, img, bounds, block_r_idx, block_c_idx,
                      intensity_range, block_cdf)
    # South-West
    bounds = ((r₁ + 1):height, 1:(c₀ - 1))
    block_r_idx = rblocks
    block_c_idx = 1
    transform_corner!(out, img, bounds, block_r_idx , block_c_idx,
                      intensity_range, block_cdf)
    # South-East
    bounds = ((r₁ + 1):height, (c₁ + 1):width)
    block_r_idx = rblocks
    block_c_idx = cblocks
    transform_corner!(out, img, bounds, block_r_idx, block_c_idx,
                      intensity_range, block_cdf)
end

function transform_interior!(out, img, bounds, block_centroids, block_dimensions, intensity_range, block_cdf)
    rows, cols = bounds
    block_centroid_r, block_centroid_c = block_centroids
    block_width, block_height = block_dimensions
    minval, maxval = intensity_range
    inv_block_height = 1 / block_height
    inv_block_width = 1 / block_width
    for r in rows
        for c in cols
            rᵢ = round(Int, r * inv_block_height)
            cᵢ = round(Int, c * inv_block_width)
            t = (r - block_centroid_r[rᵢ]) / block_height
            u = (c - block_centroid_c[cᵢ]) / block_width
            T₁ = apply_cdf_transform(img[r,c], minval, maxval, block_cdf[rᵢ, cᵢ]...)
            T₂ = apply_cdf_transform(img[r,c], minval, maxval, block_cdf[rᵢ + 1, cᵢ]...)
            T₃ = apply_cdf_transform(img[r,c], minval, maxval, block_cdf[rᵢ + 1, cᵢ + 1]...)
            T₄ = apply_cdf_transform(img[r,c], minval, maxval, block_cdf[rᵢ, cᵢ + 1]...)
            interpolated_val = (1 - t)*(1 - u)*T₁ + t*(1 - u)*T₂  + t*u*T₃ + (1 - t)*u*T₄
            out[r,c] = eltype(img) <: Integer ? ceil(interpolated_val) : interpolated_val
        end
    end
end


function transform_vertical_strip!(out, img, bounds, block_centroid_r, cᵢ, block_height, intensity_range, block_cdf)
    rows, cols = bounds
    minval, maxval = intensity_range
    inv_block_height = 1 / block_height
    for r in rows
        for c in cols
            rᵢ = round(Int, r * inv_block_height)
            t = (r - block_centroid_r[rᵢ]) / block_height
            T₁ = apply_cdf_transform(img[r,c], minval, maxval, block_cdf[rᵢ, cᵢ]...)
            T₂ = apply_cdf_transform(img[r,c], minval, maxval, block_cdf[rᵢ + 1, cᵢ]...)
            interpolated_val = (1 - t)*T₁ + t*T₂
            out[r,c] = eltype(img) <: Integer ? ceil(interpolated_val) : interpolated_val
        end
    end
end


function transform_horizontal_strip!(out, img, bounds, block_centroid_c, rᵢ, block_width, intensity_range, block_cdf)
    rows, cols = bounds
    minval, maxval = intensity_range
    inv_block_width = 1 / block_width
    for r in rows
        for c in cols
            cᵢ = round(Int, c * inv_block_width)
            u = (c - block_centroid_c[cᵢ]) / block_width
            T₁ = apply_cdf_transform(img[r,c], minval, maxval, block_cdf[rᵢ, cᵢ]...)
            T₂ = apply_cdf_transform(img[r,c], minval, maxval, block_cdf[rᵢ, cᵢ + 1]...)
            interpolated_val = (1 - u)*T₁ + u*T₂
            out[r,c] = eltype(img) <: Integer ? ceil(interpolated_val) : interpolated_val
        end
    end
end


function transform_corner!(out, img, bounds, rᵢ, cᵢ, intensity_range, block_cdf)
    rows, cols = bounds
    minval, maxval = intensity_range
    for r in rows
        for c in cols
            val = apply_cdf_transform(img[r,c], minval, maxval, block_cdf[rᵢ, cᵢ]...)
            out[r,c] = eltype(img) <: Integer ? ceil(val) : val
        end
    end
end
