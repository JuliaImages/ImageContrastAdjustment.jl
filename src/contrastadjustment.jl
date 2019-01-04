

function partition_interval(nbins::Integer, minval::Real, maxval::Real)
    edges = (range(0,step=(maxval - minval) / nbins,length=nbins)) .+ minval
end

function partition_interval(nbins::Integer, minval::AbstractGray, maxval::AbstractGray)
    partition_interval(nbins, minval.val, maxval.val)
end

"""
```
edges, count = build_histogram(img, nbins)
edges, count = build_histogram(img, nbins; minval, maxval)
edges, count = build_histogram(img, edges)
```

Generates a histogram for the image over `nbins` spread between `[minval, maxval]`.
Color images are automatically converted to grayscale.

# Output
Returns `edges` which is a `AbstractRange` type that specifies how the  interval
`[minval, maxval]` is divided into bins, and an array `count` which records the
concomitant bin frequencies. In particular, `count` has the following
properties:

* `count[0]` is the number satisfying `x < edges[1]`
* `count[i]` is the number of values `x` that satisfy `edges[i] <= x < edges[i+1]`
* `count[end]` is the number satisfying `x >= edges[end]`.
* `length(count) == length(edges)+1`.

# Details

One can consider a histogram as a piecewise-constant model of a probability
density function ``f`` [1]. Suppose that ``f`` has support on some interval ``I =
[a,b]``.  Let ``m`` be an integer and ``a = a_1 < a_2 < \\ldots < a_m < a_{m+1} =
b`` a sequence of real numbers. Construct a sequence of intervals

```math
I_1 = [a_1,a_2], I_2 = (a_2, a_3], \\ldots, I_{m} = (a_m,a_{m+1}]
```

which partition ``I`` into subsets ``I_j`` ``(j = 1, \\ldots, m)`` on which
``f`` is constant. These subsets satisfy ``I_i \\cap I_j = \\emptyset, \\forall
i \\neq j``, and are commonly referred to as *bins*. Together they encompass the
entire range of data values such that ``\\sum_j |I_j | = | I |``. Each bin has
width ``w_j = |I_j| = a_{j+1} - a_j`` and height ``h_j`` which is the constant
probability density over the region of the bin. Integrating the constant
probability density over the width of the bin ``w_j`` yields a probability mass
of ``\\pi_j = h_j w_j`` for the bin.

For a sample ``x_1, x_2, \\ldots, x_N``, let


```math
n_j = \\sum_{n = 1}^{N}\\mathbf{1}_{(I_j)}(x_n),
\\quad \\text{where} \\quad
\\mathbf{1}_{(I_j)}(x) =
\\begin{cases}
 1 & \\text{if} \\; x \\in I_j,\\\\
 0 & \\text{otherwise},
\\end{cases},
```
represent the number of samples falling into the interval ``I_j``. An estimate
for the probability mass of the ``j``th bin is given by the relative frequency
``\\hat{\\pi} = \\frac{n_j}{N}``, and the histogram estimator of the probability
density function is defined as
```math
\\begin{aligned}
\\hat{f}_n(x)  & = \\sum_{j = 1}^{m}\\frac{n_j}{Nw_j} \\mathbf{1}_{(I_j)}(x) \\\\
& = \\sum_{j = 1}^{m}\\frac{\\hat{\\pi}_j}{w_j} \\mathbf{1}_{(I_j)}(x) \\\\
& = \\sum_{j = 1}^{m}\\hat{h}_j \\mathbf{1}_{(I_j)}(x).
\\end{aligned}
```

The function ``\\hat{f}_n(x)`` is a genuine density estimator because ``\\hat{f}_n(x)  \\ge 0`` and
```math
\\begin{aligned}
\\int_{-\\infty}^{\\infty}\\hat{f}_n(x) \\operatorname{d}x & = \\sum_{j=1}^{m} \\frac{n_j}{Nw_j} w_j \\\\
& = 1.
\\end{aligned}
```

# Options
Various options for the parameters of this function are described in more detail
below.

## Choices for `nbins`
You can specify the number of discrete bins for the histogram. When specifying
the number of bins consider the maximum number of graylevels that your image
type supports. For example, with an image of type `N0f8` there is a maximum
of 256 possible graylevels. Hence, if you request more than 256 bins for
that type of image you should expect to obtain zero counts for numerous bins.


## Choices for `minval`
You have the option to specify the lower bound of the interval over which the
histogram will be computed.  If `minval` is not specified then the minimum
value present in the image is taken as the lower bound.

## Choices for `maxval`
You have the option to specify the upper bound of the interval over which the
histogram will be computed.  If `maxval` is not specified then the maximum
value present in the image is taken as the upper bound.

## Choices for `edges`
If you do not designate the number of bins, nor the lower or upper bound of the
interval, then you have the option to directly stipulate how the intervals will
be divided by specifying a `AbstractRange` type.

# Example

Compute the histogram of a grayscale image.
```julia

using TestImages, FileIO, ImageView

img =  testimage("mandril_gray");
edges, counts  = build_histogram(img, 256, minval = 0, maxval = 1)
```

Given a color image, compute the histogram of the red channel.
```julia
img = testimage("mandrill")
r = red.(img)
edges, counts  = build_histogram(r, 256, minval = 0, maxval = 1)
```

# References
[1] E. Herrholz, "Parsimonious Histograms," Ph.D. dissertation, Inst. of Math. and Comp. Sci., University of Greifswald, Greifswald, Germany, 2011.

# See also:

| Operation                     	| Function Name                                                                                                                                   	| In-place Variant                                                                                                                                  	|
|-------------------------------	|-------------------------------------------------------------------------------------------------------------------------------------------------	|---------------------------------------------------------------------------------------------------------------------------------------------------	|
| Histogram Equalization        	| [`adjust_histogram`](@ref adjust_histogram(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) 	| [`adjust_histogram!`](@ref adjust_histogram!(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) 	|
| Midway Histogram Equalization 	| [`adjust_histogram`](@ref adjust_histogram(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))                                  	| [`adjust_histogram!`](@ref adjust_histogram!(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))                                  	|
| Histogram Matching            	| [`adjust_histogram`](@ref adjust_histogram(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))                                            	| [`adjust_histogram!`](@ref adjust_histogram!(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))                                            	|
| Gamma Correction              	| [`adjust_histogram`](@ref adjust_histogram(::GammaCorrection,  ::AbstractArray, ::Real))                                                        	| [`adjust_histogram!`](@ref adjust_histogram!(::GammaCorrection,  ::AbstractArray, ::Real))                                                        	|
| Linear Stretching             	| [`adjust_histogram`](@ref adjust_histogram(::LinearStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))       	| [`adjust_histogram!`](@ref adjust_histogram!(::LinearStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))       	|
| Contrast Stretching           	| [`adjust_histogram`](@ref adjust_histogram(::ContrastStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))     	| [`adjust_histogram!`](@ref adjust_histogram!(::ContrastStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))     	|

"""
function build_histogram(img::AbstractArray, nbins::Integer = 256; minval::Union{Real,AbstractGray}, maxval::Union{Real,AbstractGray})
    edges = partition_interval(nbins, minval, maxval)
    build_histogram(img, edges)
end

function build_histogram(img::AbstractArray, nbins::Integer = 256)
    build_histogram(img, nbins, minval =  minfinite(img), maxval = maxfinite(img))
end

function build_histogram(img::AbstractArray{T}, edges::AbstractRange) where {T<:Color3}
    build_histogram(Gray.(img),edges)
end

function build_histogram(img::AbstractArray, edges::AbstractRange)
    Base.has_offset_axes(edges) && throw(ArgumentError("edges must be indexed starting with 1"))
    lb = first(axes(edges,1))-1
    ub = last(axes(edges,1))
    first_edge = first(edges)
    step_size = step(edges)
    counts = fill(0, lb:ub)
    for val in img
         if isnan(val)
             continue
         else
            if val >= edges[end]
                counts[end] += 1
                continue
            elseif val < first(edges)
                counts[lb] += 1
            else
                index = Int(Base.div(val-first_edge,step_size)) + 1
                counts[index] += 1
            end
        end
    end
    edges, counts
end

"""
```
adjust_histogram(Equalization(),img, nbins)
adjust_histogram(Equalization(),img, nbins; minval, maxval)
```

Returns a histogram equalized image with a granularity of `nbins` number of bins.

# Details

Histogram equalization was initially conceived to  improve the contrast in a
single-channel grayscale image. The method transforms the
distribution of the intensities in an image so that they are as uniform as
possible [1]. The natural justification for uniformity
is that the image has better contrast  if the intensity levels of an image span
a wide range on the intensity scale. As it turns out, the necessary
transformation is a mapping based on the cumulative histogram.

One can consider an ``L``-bit single-channel ``I \\times J`` image with gray
values in the set ``\\{0,1,\\ldots,L-1 \\}``, as a collection of independent and
identically distributed random variables. Specifically, let the sample space
``\\Omega`` be the set of all ``IJ``-tuples ``\\omega
=(\\omega_{11},\\omega_{12},\\ldots,\\omega_{1J},\\omega_{21},\\omega_{22},\\ldots,\\omega_{2J},\\omega_{I1},\\omega_{I2},\\ldots,\\omega_{IJ})``,
where each ``\\omega_{ij} \\in \\{0,1,\\ldots, L-1 \\}``. Furthermore, impose a
probability measure on ``\\Omega`` such that the functions ``\\Omega \\ni
\\omega \\to \\omega_{ij} \\in \\{0,1,\\ldots,L-1\\}`` are independent and
identically distributed.

One can then regard an image as a matrix of random variables ``\\mathbf{G} =
[G_{i,j}(\\omega)]``, where each function ``G_{i,j}: \\Omega \\to \\mathbb{R}``
is defined by
```math
G_{i,j}(\\omega) = \\frac{\\omega_{ij}}{L-1},
```
and each ``G_{i,j}`` is distributed according to some unknown density ``f_{G}``.
While ``f_{G}`` is unknown, one can approximate it with a normalized histogram
of gray levels,

```math
\\hat{f}_{G}(v)= \\frac{n_v}{IJ},
```
where
```math
n_v = \\left | \\left\\{(i,j)\\, |\\,  G_{i,j}(\\omega)  = v \\right \\} \\right |
```
represents the number of times a gray level with intensity ``v`` occurs in
``\\mathbf{G}``. To transform the distribution of the intensities so that
they are as uniform as possible one needs to find a mapping ``T(\\cdot)`` such
that ``T(G_{i,j}) \\thicksim U ``. The required mapping turns out to be the
cumulative distribution function (CDF) of the empirical density
``\\hat{f}_{G}``,
```math
 T(G_{i,j}) = \\int_0^{G_{i,j}}\\hat{f}_{G}(w)\\mathrm{d} w.
```

# Options

Various options for the parameters of this function are described in more detail
below.

## Choices for `img`

The `adjust_histogram(Equalization(),...)` function can handle a variety of
input types.  The type of the returned image matches the input type.

For colored images, the input is converted to
[YIQ](https://en.wikipedia.org/wiki/YIQ) type and the Y channel is equalized.
This is the combined with the I and Q channels and the resulting image converted
to the same type as the input.

## Choices for `nbins`

You can specify the total number of bins in the histogram.

## Choices for `minval` and `maxval`

If `minval` and `maxval` are specified then intensities are equalized to the range
[`minval`, `maxval`]. The default values are 0 and 1.

# Example

```julia

using TestImages, FileIO, ImageView

img =  testimage("mandril_gray")
imgeq = adjust_histogram(Equalization(),img,256, minval = 0, maxval = 1)

imshow(img)
imshow(imgeq)
```

# References
1. R. C. Gonzalez and R. E. Woods. *Digital Image Processing (3rd Edition)*.  Upper Saddle River, NJ, USA: Prentice-Hall,  2006.

# See also:

| Operation                     	| Function Name                                                                                                                                   	| In-place Variant                                                                                                                                  	|
|-------------------------------	|-------------------------------------------------------------------------------------------------------------------------------------------------	|---------------------------------------------------------------------------------------------------------------------------------------------------	|
| Histogram Construction        	| [`build_histogram`](@ref build_histogram(::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))                   	|                                                                                                                                                   	|
| Midway Histogram Equalization 	| [`adjust_histogram`](@ref adjust_histogram(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))                                  	| [`adjust_histogram!`](@ref adjust_histogram!(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))                                  	|
| Histogram Matching            	| [`adjust_histogram`](@ref adjust_histogram(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))                                            	| [`adjust_histogram!`](@ref adjust_histogram!(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))                                            	|
| Gamma Correction              	| [`adjust_histogram`](@ref adjust_histogram(::GammaCorrection,  ::AbstractArray, ::Real))                                                        	| [`adjust_histogram!`](@ref adjust_histogram!(::GammaCorrection,  ::AbstractArray, ::Real))                                                        	|
| Linear Stretching             	| [`adjust_histogram`](@ref adjust_histogram(::LinearStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))       	| [`adjust_histogram!`](@ref adjust_histogram!(::LinearStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))       	|
| Contrast Stretching           	| [`adjust_histogram`](@ref adjust_histogram(::ContrastStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))     	| [`adjust_histogram!`](@ref adjust_histogram!(::ContrastStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))     	|

"""
function adjust_histogram(operation::Equalization, img::AbstractArray, nbins::Integer = 256; minval::Union{Real,AbstractGray} = 0, maxval::Union{Real,AbstractGray} = 1)
    adjust_histogram!(Equalization(), copy(img), nbins, minval = minval, maxval = maxval)
end


function adjust_histogram(operation::Equalization, img::AbstractArray{T}, nbins::Integer = 256; minval::Union{Real,AbstractGray} = 0, maxval::Union{Real,AbstractGray} = 1) where {T<:Color3}
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    adjust_histogram!(Equalization(), view(yiq_view,1,:,:), nbins, minval = minval, maxval = maxval)
    convert.(T, yiq)
end


"""
```
adjust_histogram!(Equalization(),img, nbins)
adjust_histogram!(Equalization(),img, nbins; minval, maxval)
```

Same as [`adjust_histogram`](@ref adjust_histogram(::Equalization, ::AbstractArray, ::Integer, ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) except that it modifies the image that was passed as an argument.
"""
function adjust_histogram!(operation::Equalization, img::AbstractArray, nbins::Integer; minval::Union{Real,AbstractGray} = 0, maxval::Union{Real,AbstractGray} = 1)
    edges, histogram = build_histogram(img, nbins, minval = minval, maxval = maxval)
    lb = first(axes(histogram,1))
    ub = last(axes(histogram,1))
    N = length(img)
    cdf = cumsum(histogram[lb:ub]/N)
    transform_density!(img, edges, cdf, minval, maxval)
end

function transform_density!(img::AbstractArray,edges::AbstractArray, cdf::AbstractArray, minval::Union{Real,AbstractGray}, maxval::Union{Real,AbstractGray})
    first_edge = first(edges)
    step_size = step(edges)
    T = eltype(img)
    map!(img,img) do val
        if val >= edges[end]
            newval = cdf[end]
        elseif val < first_edge
            newval = first(cdf)
        else
            index = Int(Base.div(val-first_edge,step_size)) + 1
            newval = cdf[index]
        end
        # Scale the new intensity value to so that it lies in the range [minval, maxval].
        if T <: Integer
            newval = ceil(minval + ((newval - first(cdf)) * (maxval - minval) / (cdf[end] - first(cdf))))
        else
            newval = minval + ((newval - first(cdf)) * (maxval - minval) / (cdf[end] - first(cdf)))
        end
    end
end


"""
```
adjust_histogram(Matching(),img, targetimg, nbins)
adjust_histogram(Matching(),img, targetimg, edges)
```

Returns a histogram matched image with a granularity of `nbins` number of bins.
The first argument `img` is the image to be matched, and the second argument
`targetimg` is the image having the desired histogram to be matched to.

# Details
The purpose of histogram matching is to transform the intensities in a source
image so that the intensities distribute according to the histogram of a
specified target image. If one interprets histograms as piecewise-constant
models of probability density functions (see [`build_histogram`](@ref
build_histogram(::AbstractArray, ::Integer, ::Union{Real,AbstractGray},
::Union{Real,AbstractGray}))), then the histogram matching task can be modelled
as the problem of transforming one probability distribution into another [1].
It turns out that the solution to this transformation problem involves the
cumulative and inverse cumulative distribution functions of the source and
target probability density functions.

In particular, let the random variables ``x \\thicksim p_{x} `` and ``z
\\thicksim p_{z}``  represent an intensity in the source and target image
respectively, and let

```math
 S(x) = \\int_0^{x}p_{x}(w)\\mathrm{d} w \\quad \\text{and} \\quad
 T(z) = \\int_0^{z}p_{z}(w)\\mathrm{d} w
```
represent their concomitant cumulative distribution functions. Then the
sought-after mapping ``Q(\\cdot)`` such that ``Q(x) \\thicksim p_{z} `` is given
by

```math
Q(x) =  T^{-1}\\left( S(x) \\right),
```

where ``T^{-1}(y) = \\operatorname{min} \\{ x \\in \\mathbb{R} : y \\leq T(x)
\\}`` is the inverse cumulative distribution function of ``T(x)``.

The mapping suggests that one can conceptualize histogram matching as performing
histogram equalization on the source and target image and relating the two
equalized histograms. Refer to [`adjust_histogram`](@ref
adjust_histogram(::Equalization, ::AbstractArray, ::Integer,
::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) for more details on
histogram equalization.

# Options

Various options for the parameters of this function are described in more detail
below.

## Choices for `img` and `targetimg`

The `adjust_histogram(Matching(),...)` function can handle a variety of input
types. The type of the returned image matches the input type.

For colored images, the inputs are converted to
[YIQ](https://en.wikipedia.org/wiki/YIQ)  type and the distributions of the Y
channels are matched. The modified Y channel is then combined with the I and Q
channels and the resulting image converted to the same type as the input.

## Choices for `nbins`

You can specify the total number of bins in the histogram. If you do not
specify the number of bins then a default value of 256 bins is utilized.

## Choices for `edges`

If you do not designate the number of bins, then you have the option to directly
stipulate how the intervals will be divided by specifying a `AbstractRange`
type.

# Example

```julia
using Images, TestImages, ImageView

img_source = testimage("mandril_gray")
img_target = adjust_gamma(img_source, 1/2)
img_transformed = adjust_histogram(Matching(),img_source, img_target)
#=
    A visual inspection confirms that img_transformed resembles img_target
    much more closely than img_source.
=#
imshow(img_source)
imshow(img_target)
imshow(img_transformed)
```

# References
1. W. Burger and M. J. Burge. *Digital Image Processing*. Texts in Computer Science, 2016. [doi:10.1007/978-1-4471-6684-9](https://doi.org/10.1007/978-1-4471-6684-9)


# See also:

| Operation                     	| Function Name                                                                                                                                   	| In-place Variant                                                                                                                                  	|
|-------------------------------	|-------------------------------------------------------------------------------------------------------------------------------------------------	|---------------------------------------------------------------------------------------------------------------------------------------------------	|
| Histogram Construction        	| [`build_histogram`](@ref build_histogram(::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))                   	|                                                                                                                                                   	|
| Histogram Equalization        	| [`adjust_histogram`](@ref adjust_histogram(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) 	| [`adjust_histogram!`](@ref adjust_histogram!(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) 	|
| Midway Histogram Equalization 	| [`adjust_histogram`](@ref adjust_histogram(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))                                  	| [`adjust_histogram!`](@ref adjust_histogram!(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))                                  	|
| Gamma Correction              	| [`adjust_histogram`](@ref adjust_histogram(::GammaCorrection,  ::AbstractArray, ::Real))                                                        	| [`adjust_histogram!`](@ref adjust_histogram!(::GammaCorrection,  ::AbstractArray, ::Real))                                                        	|
| Linear Stretching             	| [`adjust_histogram`](@ref adjust_histogram(::LinearStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))       	| [`adjust_histogram!`](@ref adjust_histogram!(::LinearStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))       	|
| Contrast Stretching           	| [`adjust_histogram`](@ref adjust_histogram(::ContrastStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))     	| [`adjust_histogram!`](@ref adjust_histogram!(::ContrastStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))     	|

"""
function adjust_histogram(operation::Matching, img::AbstractArray, targetimg::AbstractArray, nbins::Integer = 256)
    adjust_histogram!(Matching(), copy(img), targetimg, nbins)
end


function adjust_histogram(operation::Matching, img::AbstractArray, targetimg::AbstractArray, edges::AbstractRange)
    adjust_histogram!(Matching(), copy(img), targetimg, edges)
end


function adjust_histogram!(operation::Matching, img::AbstractArray{T}, targetimg::AbstractArray{T}, edges::AbstractRange ) where T <: Color3
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    adjust_histogram!(Matching(),view(yiq_view,1,:,:), map(c->YIQ(c).y, targetimg), edges)
    convert.(T, yiq)
end


function adjust_histogram!(operation::Matching, img::AbstractArray{T}, targetimg::AbstractArray{T}, nbins::Integer = 256 ) where T <: Color3
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    adjust_histogram!(Matching(),view(yiq_view,1,:,:), map(c->YIQ(c).y, targetimg), nbins)
    convert.(T, yiq)
end


function adjust_histogram!(operation::Matching, img::AbstractArray, targetimg::AbstractArray, edges::AbstractRange )
    edges, pdf, target_pdf = construct_pdfs(img, targetimg, edges)
    match_pdf!(Matching(), img, edges, pdf, target_pdf)
end

"""
```
adjust_histogram!(Matching(),img, targetimg, nbins)
adjust_histogram!(Matching(),img, targetimg, edges)
```

Same as  [`adjust_histogram`](@ref adjust_histogram(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))  except that it modifies the image that was passed as an argument.
"""
function adjust_histogram!(operation::Matching, img::AbstractArray, targetimg::AbstractArray, nbins::Integer = 256 )
    edges, pdf, target_pdf = construct_pdfs(img, targetimg, nbins)
    match_pdf!(Matching(), img, edges, pdf, target_pdf)
end

function construct_pdfs(img::AbstractArray, targetimg::AbstractArray, edges::AbstractRange)
    _, histogram = build_histogram(img, edges)
    _, target_histogram = build_histogram(targetimg, edges)
    return edges, histogram / sum(histogram), target_histogram / sum(target_histogram)
end

function construct_pdfs(img::AbstractArray, targetimg::AbstractArray, nbins::Integer = 256)
    if eltype(img) <: AbstractGray
        imin, imax = 0, 1
    else
        imin, imax = min(minfinite(img), minfinite(targetimg)), max(maxfinite(img), maxfinite(targetimg))
    end
    edges, histogram = build_histogram(img, nbins, minval = imin, maxval = imax)
    _, target_histogram = build_histogram(targetimg, edges)
    return edges, histogram / sum(histogram), target_histogram / sum(target_histogram)
end

function lookup_icdf(cdf::AbstractArray, targetcdf::AbstractArray)
    lookup_table = zeros(Int, length(cdf))
    i = 1
    for j = 1:length(cdf)
        p = cdf[j]
        while i < length(targetcdf) && targetcdf[i+1] <= p
            i += 1
        end
        lookup_table[j] = i
    end
    lookup_table
end

function match_pdf!(operation::Matching, img::AbstractArray, edges::AbstractArray, pdf::AbstractArray, target_pdf::AbstractArray)
    cdf = parent(cumsum(pdf))
    target_cdf = parent(cumsum(target_pdf))
    # Precompute the inverse cummulative distribution function of target_cdf.
    lookup_table = lookup_icdf(cdf, target_cdf)
    # Transform the intensities in img so that they are distributed according
    # to the distribution of the target_histogram.
    T = eltype(img)
    step_size = step(edges)
    first_edge = first(edges)
    last_edge = last(edges)
    map!(img,img) do val
        if isnan(val)
            return val
        else
            if val >= last_edge
                newval = edges[last(lookup_table)-1] + step_size
            elseif val < first_edge
                newval = edges[first(lookup_table)]
            else
                index = Int(Base.div(val-first_edge,step_size)) + 1
                newval = edges[lookup_table[index]]
            end
            return T <: Integer ? ceil(newval) : newval
        end
    end
end

"""
```
adjust_histogram(GammaCorrection(),img, gamma)
```

Returns a gamma corrected image.

# Details


Gamma correction is a non-linear  transformation given by the relation
```math
f(x) = x^\\gamma \\quad \\text{for} \\; x \\in \\mathbb{R}, \\gamma > 0.
```
It is called a *power law* transformation because one quantity varies as a power
of another quantity.

Gamma correction has historically been used to preprocess
an image to compensate for the fact that the intensity of light generated by a
physical device is not usually a linear function of the applied signal but
instead follows a power law [1]. For example, for many Cathode Ray Tubes (CRTs) the
emitted light intensity on the display is approximately equal to the voltage
raised to the power of γ, where γ ∈ [1.8, 2.8]. Hence preprocessing a raw image with
an exponent of 1/γ  would have ensured a linear response to brightness.

Research in psychophysics has also established an [empirical  power law
](https://en.wikipedia.org/wiki/Stevens%27s_power_law)  between light intensity and perceptual
brightness. Hence, gamma correction often serves as a useful image enhancement
tool.


# Options

Various options for the parameters of this function are described in more detail
below.

## Choices for `img`

The function can handle a variety of input types. The returned
image depends on the input type.

For colored images, the input is converted to YIQ type and the Y channel is
gamma corrected. This is the combined with the I and Q channels and the
resulting image converted to the same type as the input.

## Choice for `gamma`

The `gamma` value must be a non-zero positive number. A `gamma` value less than
one will yield a brighter image whereas a value greater than one will produce a
darker image. If left unspecified a default value of one is assumed.

# Example

```julia
using ImageContrastAdjustment, ImageView

# Create an example image consisting of a linear ramp of intensities.
n = 32
intensities = 0.0:(1.0/n):1.0
img = repeat(intensities, inner=(20,20))'

# Brighten the dark tones.
imgadj = adjust_histogram(GammaCorrection(), img, 1/2)

# Display the original and adjusted image.
imshow(img)
imshow(imgadj)
```

# References
1. W. Burger and M. J. Burge. *Digital Image Processing*. Texts in Computer Science, 2016. [doi:10.1007/978-1-4471-6684-9](https://doi.org/10.1007/978-1-4471-6684-9)

# See also:

| Operation                     	| Function Name                                                                                                                                   	| In-place Variant                                                                                                                                  	|
|-------------------------------	|-------------------------------------------------------------------------------------------------------------------------------------------------	|---------------------------------------------------------------------------------------------------------------------------------------------------	|
| Histogram Construction        	| [`build_histogram`](@ref build_histogram(::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))                   	|                                                                                                                                                   	|
| Histogram Equalization        	| [`adjust_histogram`](@ref adjust_histogram(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) 	| [`adjust_histogram!`](@ref adjust_histogram!(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) 	|
| Midway Histogram Equalization 	| [`adjust_histogram`](@ref adjust_histogram(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))                                  	| [`adjust_histogram!`](@ref adjust_histogram!(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))                                  	|
| Histogram Matching            	| [`adjust_histogram`](@ref adjust_histogram(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))                                            	| [`adjust_histogram!`](@ref adjust_histogram!(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))                                            	|
| Linear Stretching             	| [`adjust_histogram`](@ref adjust_histogram(::LinearStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))       	| [`adjust_histogram!`](@ref adjust_histogram!(::LinearStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))       	|
| Contrast Stretching           	| [`adjust_histogram`](@ref adjust_histogram(::ContrastStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))     	| [`adjust_histogram!`](@ref adjust_histogram!(::ContrastStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))     	|

"""
function adjust_histogram(operation::GammaCorrection, img::AbstractArray, gamma::Real= 1.0)
    adjust_histogram!(GammaCorrection(), copy(img), gamma)
end


function adjust_histogram(operation::GammaCorrection, img::AbstractArray{T}, gamma::Real = 1.0) where {T<:Color3}
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    adjust_histogram!(GammaCorrection(), view(yiq_view,1,:,:), gamma)
    convert.(T, yiq)
end


"""
```
adjust_histogram!(GammaCorrection(),img, gamma)
```

Same as [`adjust_histogram`](@ref adjust_histogram(::GammaCorrection, ::AbstractArray, ::Real)) except that it modifies the image that was passed as an argument.
"""
function adjust_histogram!(operation::GammaCorrection, img::AbstractArray, gamma::Real = 1.0)
    γ = Float64(gamma)
    T = eltype(img)
    map!(img,img) do val
        if isnan(val)
            return val
        else
            return  T <: Integer ? round(Int, val^γ) : val^γ
        end
    end
end

function adjust_histogram!(operation::GammaCorrection, img::AbstractArray{Gray{T}}, gamma::Real = 1.0) where {T<:FixedPointNumbers.Normed}
    γ = Float64(gamma)
    # Create a lookup-table for the gamma transformation of the grayvalues.
    raw_type = FixedPointNumbers.rawtype(T)
    table = zeros(T, typemax(raw_type) + 1)
    for i in zero(raw_type):typemax(raw_type)
        table[i + 1] = T((i / typemax(raw_type))^γ)
    end
    # Map the pixels to their new grayvalues.
    map!(img,img) do p
        table[p.val.i + 1]
    end
end


"""
```
adjust_histogram(LinearStretching(), img; minval = 0, maxval = 1)
```

Returns an image where the range of the intensities spans the interval [`minval`, `maxval`].

# Details


Linear stretching (also called *normalization*) is a contrast enhancing
transformation that is used to modify the dynamic range of the image. In
particular, suppose that the input image has gray values in the range [A,B] and
one wishes to change the dynamic range to [a,b] using a linear mapping, then the
necessary transformation is given by the relation

```math
f(x) = (x-A) \\frac{b-a}{B-A} + a.
```

# Options

Various options for the parameters of this function are described in more detail
below.

## Choices for `img`

The function can handle a variety of input types. The returned
image depends on the input type.

For colored images, the input is converted to the
[YIQ](https://en.wikipedia.org/wiki/YIQ)  type and the intensities of the Y
channel are stretched to the specified range. The modified Y channel is then
combined with the I and Q channels and the resulting image converted to the same
type as the input.

## Choices for `minval` and `maxval`

If minval and maxval are specified then intensities are mapped to the range
[`minval`, `maxval`]. The default values are 0 and 1.

# Example

```julia
using ImageContrastAdjustment, ImageView, TestImages, Images

img = testimage("mandril_gray")
imgo = adjust_histogram(LinearStretching(),img, minval = 0, maxval = 1)

```

# References
1. W. Burger and M. J. Burge. *Digital Image Processing*. Texts in Computer Science, 2016. [doi:10.1007/978-1-4471-6684-9](https://doi.org/10.1007/978-1-4471-6684-9)

# See also:

| Operation                     	| Function Name                                                                                                                                   	| In-place Variant                                                                                                                                  	|
|-------------------------------	|-------------------------------------------------------------------------------------------------------------------------------------------------	|---------------------------------------------------------------------------------------------------------------------------------------------------	|
| Histogram Construction        	| [`build_histogram`](@ref build_histogram(::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))                   	|                                                                                                                                                   	|
| Histogram Equalization        	| [`adjust_histogram`](@ref adjust_histogram(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) 	| [`adjust_histogram!`](@ref adjust_histogram!(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) 	|
| Midway Histogram Equalization 	| [`adjust_histogram`](@ref adjust_histogram(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))                                  	| [`adjust_histogram!`](@ref adjust_histogram!(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))                                  	|
| Histogram Matching            	| [`adjust_histogram`](@ref adjust_histogram(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))                                            	| [`adjust_histogram!`](@ref adjust_histogram!(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))                                            	|
| Gamma Correction              	| [`adjust_histogram`](@ref adjust_histogram(::GammaCorrection,  ::AbstractArray, ::Real))                                                        	| [`adjust_histogram!`](@ref adjust_histogram!(::GammaCorrection,  ::AbstractArray, ::Real))                                                        	|
| Contrast Stretching           	| [`adjust_histogram`](@ref adjust_histogram(::ContrastStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))     	| [`adjust_histogram!`](@ref adjust_histogram!(::ContrastStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))     	|
"""
function adjust_histogram(operation::LinearStretching, img::AbstractArray; minval::Union{Real,AbstractGray} = 0.0, maxval::Union{Real,AbstractGray} = 1.0)
    adjust_histogram!(LinearStretching(), copy(img), minval = minval, maxval = maxval)
end

"""
```
adjust_histogram!(LinearStretching(), img; minval = 0.0, maxval = 1.0)
```

Same as [`adjust_histogram`](@ref adjust_histogram(::LinearStretching, ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) except that it modifies the image that was passed as an argument.
"""
function adjust_histogram!(operation::LinearStretching, img::AbstractArray; minval::Union{Real,AbstractGray} = 0.0, maxval::Union{Real,AbstractGray} = 1.0)
    src_minval = minfinite(img)
    src_maxval = maxfinite(img)
    T = eltype(img)
    map!(img,img) do val
        if isnan(val)
            return val
        else
            newval = linear_stretch(val, src_minval, src_maxval, minval, maxval)
            return  T <: Integer ? round(Int, newval ) : newval
        end
    end
end

function adjust_histogram(operation::LinearStretching, img::AbstractArray{T}; minval::Union{Real,AbstractGray} = 0.0, maxval::Union{Real,AbstractGray} = 1.0) where {T<:Color3}
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    adjust_histogram!(LinearStretching(),view(yiq_view,1,:,:), minval = minval, maxval = maxval)
    convert.(T, yiq)
end

function linear_stretch(x, A, B, a, b)
    return (x-A) * ((b-a)/(B-A)) + a
end



"""
```
adjust_histogram(ContrastStretching(), img; t = 0.5,  slope = 1.0)
```

Returns an image where intensities intensities below `t` are compressed
into a narrower range of dark intensities, and values above `t` are compressed
into a narrower band of light intensities.

# Details

Contrast stretching is a transformation that  enhances or reduces (for `slope` >
1 or < 1, respectively) the contrast near saturation (0 and 1).
It is given by the relation
```math
f(x) = \\frac{1}{1 + \\left(\\frac{t}{x} \\right)^s}, \\; s \\in \\mathbb{R},
```
where ``s`` represents the `slope` argument.

# Options

Various options for the parameters of this function are described in more detail
below.

## Choices for `img`

The function can handle a variety of input types. The returned
image depends on the input type.

For colored images, the input is converted to the
[YIQ](https://en.wikipedia.org/wiki/YIQ)  type and the intensities of the Y
channel are stretched to the specified range. The modified Y channel is then
combined with the I and Q channels and the resulting image converted to the same
type as the input.

## Choice for `t`

The value of `t` needs to be in the unit interval. If left unspecified a
default value of 0.5 is utilized.

## Choice for `slope`

The value of `slope` can be any real number. If left unspecified a
default value of 1.0 is utilized.

# Example

```julia
using ImageContrastAdjustment, ImageView, Images, TestImages

img = testimage("mandril_gray")
ret = adjust_histogram(ContrastStretching(),img, t = 0.6, slope = 3)

```

# References
1. Gonzalez, R. C., Woods, R. E., & Eddins, S. L. (2004). *Digital image processing using MATLAB* (Vol. 624). Upper Saddle River, New Jersey: Pearson-Prentice-Hall.

# See also:

| Operation                     	| Function Name                                                                                                                                   	| In-place Variant                                                                                                                                  	|
|-------------------------------	|-------------------------------------------------------------------------------------------------------------------------------------------------	|---------------------------------------------------------------------------------------------------------------------------------------------------	|
| Histogram Construction        	| [`build_histogram`](@ref build_histogram(::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))                   	|                                                                                                                                                   	|
| Histogram Equalization        	| [`adjust_histogram`](@ref adjust_histogram(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) 	| [`adjust_histogram!`](@ref adjust_histogram!(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) 	|
| Midway Histogram Equalization 	| [`adjust_histogram`](@ref adjust_histogram(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))                                  	| [`adjust_histogram!`](@ref adjust_histogram!(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))                                  	|
| Histogram Matching            	| [`adjust_histogram`](@ref adjust_histogram(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))                                            	| [`adjust_histogram!`](@ref adjust_histogram!(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))                                            	|
| Gamma Correction              	| [`adjust_histogram`](@ref adjust_histogram(::GammaCorrection,  ::AbstractArray, ::Real))                                                        	| [`adjust_histogram!`](@ref adjust_histogram!(::GammaCorrection,  ::AbstractArray, ::Real))                                                        	|
| Linear Stretching             	| [`adjust_histogram`](@ref adjust_histogram(::LinearStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))       	| [`adjust_histogram!`](@ref adjust_histogram!(::LinearStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))       	|
"""
function adjust_histogram(operation::ContrastStretching, img::AbstractArray; t::Union{Real,AbstractGray} = 0.5, slope::Union{Real,AbstractGray} = 1.0)
    adjust_histogram!(ContrastStretching(), copy(img), t = t, slope = slope)
end


function adjust_histogram(operation::ContrastStretching, img::AbstractArray{T}; t::Union{Real,AbstractGray} = 0.5, slope::Union{Real,AbstractGray} = 1.0) where {T<:Color3}
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    adjust_histogram!(ContrastStretching(),view(yiq_view,1,:,:), t = t, slope = slope)
    convert.(T, yiq)
end

"""
```
adjust_histogram!(ContrastStretching(),img; t = 0.5, slope = 1.0)
```

Same as [`adjust_histogram`](@ref adjust_histogram(::ContrastStretching, ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) except that it modifies the image that was passed as an argument.
"""
function adjust_histogram!(operation::ContrastStretching, img::AbstractArray; t::Union{Real,AbstractGray} = 0.5, slope::Union{Real,AbstractGray} = 1.0)
    T = eltype(img)
    ϵ = eps(T)
    map!(img,img) do val
        if isnan(val)
            return val
        else
            newval = contrast_stretch(val, t, slope, ϵ)
            return  T <: Integer ? round(Int, newval ) : newval
        end
    end
end

function contrast_stretch(x, t, s, ϵ)
    return 1 / (1 + (t / (x+ϵ))^s)
end


"""
```
adjust_histogram(MidwayEqualization(),img1, img2, nbins)
adjust_histogram(MidwayEqualization(),img1, img2, edges)
```

Gives a pair of images the same histogram whilst maintaining as
much as possible their previous grey level dynamics.

# Details

The purpose of midway histogram equalization is to transform the intensities in
a pair of images so that the intensities distribute according to a common
"midway" distribution. The histogram representing the common distribution is
chosen so that the original  gray level dynamics of the images are preserved as
much as possible. If one interprets histograms as piecewise-constant models of
probability density functions (see [`build_histogram`](@ref
build_histogram(::AbstractArray, ::Integer, ::Union{Real,AbstractGray},
::Union{Real,AbstractGray}))), then the midway histogram equalization task can
be modeled as the problem of transforming one probability distribution into
another (see [`adjust_histogram`](@ref adjust_histogram(::Matching,::AbstractArray, ::AbstractArray, ::Integer))).
It turns out that the solution to this transformation problem involves the
cumulative and inverse cumulative distribution functions of the source and
"midway" probability density functions. In particular, let the random variables ``X_i \\thicksim p_{x_i} \\; (i = 1,2)``,
and ``Z \\thicksim p_{z}``  represent an intensity in the first, second and
"midway" image respectively, and let

```math
 S_{X_i}(x) = \\int_0^{x}p_{x_i}(w)\\mathrm{d} w \\; \\quad \\text{and} \\quad
 T_{Z}(x) = \\frac{2}{\\frac{1}{S_{X_1}(x)} + \\frac{1}{S_{X_2}(x)}}
```
represent the cumulative distribution functions of the two input images, and
their *harmonic mean*, respectively. Then the sought-after mapping
``Q_{X_i}(\\cdot)`` ``(i = 1,2)`` such that ``Q_{X_i}(x) \\thicksim p_{z} `` is
given by

```math
Q_{X_i}(x) =  T_{Z}^{-1}\\left( S_{X_i}(x) \\right),
```

where ``T_{Z}^{-1}(y) = \\operatorname{min} \\{ x \\in \\mathbb{R} : y \\leq
T_{Z}(x) \\}`` is the inverse cumulative distribution function of ``T_{Z}(x)``.


# Options

Various options for the parameters of this function are described in more detail
below.

## Choices for `img1` and `img2`

The  function can handle a variety of input
types. The type of the returned image matches the input type.

For colored images, the inputs are converted to
[YIQ](https://en.wikipedia.org/wiki/YIQ)  type and the distributions of the Y
channels are transformed according to a "midway" distribution. The modified Y
channel is then combined with the I and Q channels and the resulting image
converted to the same type as the input.

## Choices for `nbins`

You can specify the total number of bins in the histogram. If you do not
specify the number of bins then a default value of 256 bins is utilized.

## Choices for `edges`

If you do not designate the number of bins, then you have the option to directly
stipulate how the intervals will be divided by specifying a `AbstractRange`
type.

# Example

```julia
using Images, TestImages, ImageView, ImageContrastAdjustment

img = testimage("mandril_gray")

# The same image but with different intensitiy distributions
img1 = adjust_histogram(GammaCorrection(),img, 2)
img2 = adjust_histogram(GammaCorrection(),img, 1.2)

# Midway histogram equalization will transform these two images so that their
# intensity distributions are almost identical.
img1o, img2o = adjust_histogram(MidwayEqualization(),img1, img2, 256)

```

# References
1. T. Guillemot and J. Delon, “*Implementation of the Midway Image Equalization*,” Image Processing On Line, vol. 5, pp. 114–129, Jun. 2016. [doi:10.5201/ipol.2016.140](https://doi.org/10.5201/ipol.2016.140)


# See also:

| Operation                     	| Function Name                                                                                                                                   	| In-place Variant                                                                                                                                  	|
|-------------------------------	|-------------------------------------------------------------------------------------------------------------------------------------------------	|---------------------------------------------------------------------------------------------------------------------------------------------------	|
| Histogram Construction        	| [`build_histogram`](@ref build_histogram(::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))                   	|                                                                                                                                                   	|
| Histogram Equalization        	| [`adjust_histogram`](@ref adjust_histogram(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) 	| [`adjust_histogram!`](@ref adjust_histogram!(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})) 	|
| Histogram Matching            	| [`adjust_histogram`](@ref adjust_histogram(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))                                            	| [`adjust_histogram!`](@ref adjust_histogram!(::Matching, ::AbstractArray, ::AbstractArray, ::Integer))                                            	|
| Gamma Correction              	| [`adjust_histogram`](@ref adjust_histogram(::GammaCorrection,  ::AbstractArray, ::Real))                                                        	| [`adjust_histogram!`](@ref adjust_histogram!(::GammaCorrection,  ::AbstractArray, ::Real))                                                        	|
| Linear Stretching             	| [`adjust_histogram`](@ref adjust_histogram(::LinearStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))       	| [`adjust_histogram!`](@ref adjust_histogram!(::LinearStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))       	|
| Contrast Stretching           	| [`adjust_histogram`](@ref adjust_histogram(::ContrastStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))     	| [`adjust_histogram!`](@ref adjust_histogram!(::ContrastStretching,  ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray}))     	|

"""
function adjust_histogram(operation::MidwayEqualization, img1::AbstractArray, img2::AbstractArray, nbins::Integer = 256)
    adjust_histogram!(MidwayEqualization(), copy(img1), copy(img2), nbins)
end


function adjust_histogram(operation::MidwayEqualization, img1::AbstractArray, img2::AbstractArray, edges::AbstractRange)
    adjust_histogram!(MidwayEqualization(), copy(img1), copy(img2), edges)
end


function adjust_histogram!(operation::MidwayEqualization, img1::AbstractArray{T}, img2::AbstractArray{T}, edges::AbstractRange) where T <: Color3
    yiq1 = convert.(YIQ, img1)
    yiq1_view = channelview(yiq1)

    yiq2 = convert.(YIQ, img2)
    yiq2_view = channelview(yiq2)

    adjust_histogram!(MidwayEqualization(),view(yiq1_view,1,:,:), view(yiq2_view,1,:,:), edges)
    convert.(T, yiq1), convert.(T, yiq2)
end


function adjust_histogram!(operation::MidwayEqualization, img1::AbstractArray{T}, img2::AbstractArray{T}, nbins::Integer = 256) where T <: Color3
    yiq1 = convert.(YIQ, img1)
    yiq1_view = channelview(yiq1)

    yiq2 = convert.(YIQ, img2)
    yiq2_view = channelview(yiq2)

    adjust_histogram!(MidwayEqualization(),view(yiq1_view,1,:,:), view(yiq2_view,1,:,:), nbins)
    convert.(T, yiq1), convert.(T, yiq2)
end


function adjust_histogram!(operation::MidwayEqualization, img1::AbstractArray, img2::AbstractArray, edges::AbstractRange)
    edges, pdf1, pdf2 = construct_pdfs(img1, img2, edges)
    midway_pdf = zero(pdf1)
    cdf1 = cumsum(pdf1)
    cdf2 = cumsum(pdf2)
    # midway_cdf is the harmonic mean between cdf1 and cdf2.
    midway_cdf = similar(cdf1)
    for i in eachindex(cdf1)
        if cdf1[i] == 0 || cdf2[i] == 0
            midway_cdf[i] = 0
        else
            midway_cdf[i] = 2 / (1/cdf1[i] + 1/cdf2[i])
        end
    end
    cdf2pdf!(midway_pdf, midway_cdf)
    match_pdf!(Matching(), img1, edges, pdf1, midway_pdf),  match_pdf!(Matching(), img2, edges, pdf2, midway_pdf)
end

"""
```
adjust_histogram!(MidwayEqualization(),img1, img2, nbins)
adjust_histogram!(MidwayEqualization(),img1, img2, edges)
```

Same as  [`adjust_histogram`](@ref adjust_histogram(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer))  except that it modifies the images that were passed as arguments.
"""
function adjust_histogram!(operation::MidwayEqualization, img1::AbstractArray, img2::AbstractArray, nbins::Integer = 256)
    edges, pdf1, pdf2 = construct_pdfs(img1, img2, nbins)
    midway_pdf = zero(pdf1)
    cdf1 = cumsum(pdf1)
    cdf2 = cumsum(pdf2)
    # midway_cdf is the harmonic mean between cdf1 and cdf2.
    midway_cdf =  similar(cdf1)
    for i in eachindex(cdf1)
        if cdf1[i] == 0 || cdf2[i] == 0
            midway_cdf[i] = 0
        else
            midway_cdf[i] = 2 / (1/cdf1[i] + 1/cdf2[i])
        end
    end
    cdf2pdf!(midway_pdf, midway_cdf)
    match_pdf!(Matching(), img1, edges, pdf1, midway_pdf),  match_pdf!(Matching(), img2, edges, pdf2, midway_pdf)
end

function cdf2pdf!(pdf::AbstractArray, cdf::AbstractArray)
    pdf[1] = cdf[1]
    for i = length(cdf)-1:-1:2
        pdf[i] = cdf[i] - cdf[i-1]
    end
end
