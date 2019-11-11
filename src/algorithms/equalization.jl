"""
```
    Equalization <: AbstractHistogramAdjustmentAlgorithm
    Equalization(; nbins = 256, minval = 0, maxval = 1)

    adjust_histogram([T,] img, f::Equalization)
    adjust_histogram!([out,] img, f::Equalization)
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

Various options for the parameters of the `adjust_histogram` function and
`Equalization` type are described in more detail below.

## Choices for `img`

The `adjust_histogram` function can handle a variety of
input types.  By default type of the returned image matches the input type.

For colored images, the input is converted to
[YIQ](https://en.wikipedia.org/wiki/YIQ) type and the Y channel is equalized.
This is the combined with the I and Q channels and the resulting image converted
to the same type as the input.

## Choices for `nbins` in `Equalization`

You can specify the total number of bins in the histogram.

## Choices for `minval` and `maxval` in `Equalization`

If `minval` and `maxval` are specified then intensities are equalized to the range
[`minval`, `maxval`]. The default values are 0 and 1.

# Example

```julia

using TestImages, FileIO, ImageView

img =  testimage("mandril_gray")
imgeq = adjust_histogram(img, Equalization(nbins = 256, minval = 0, maxval = 1))

imshow(img)
imshow(imgeq)
```

# References
1. R. C. Gonzalez and R. E. Woods. *Digital Image Processing (3rd Edition)*.  Upper Saddle River, NJ, USA: Prentice-Hall,  2006.
"""
Base.@kwdef struct Equalization{T₁ <: Union{Real,AbstractGray},
                                T₂ <: Union{Real,AbstractGray}} <: AbstractHistogramAdjustmentAlgorithm
    nbins::Int = 256
    minval::T₁ = 0.0
    maxval::T₂ = 1.0
end

function (f::Equalization)(out::GenericGrayImage, img::GenericGrayImage)
    edges, histogram = build_histogram(img, f.nbins, minval = f.minval, maxval = f.maxval)
    lb = first(axes(histogram,1))
    ub = last(axes(histogram,1))
    N = length(img)
    cdf = cumsum(histogram[lb:ub]/N)
    out .= img
    transform_density!(out, edges, cdf, f.minval, f.maxval)
end

function (f::Equalization)(out::AbstractArray{<:Color3}, img::AbstractArray{<:Color3})
    T = eltype(img)
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    adjust_histogram!(view(yiq_view,1,:,:), f)
    out .= convert.(T, yiq)
end

(f::Equalization)(out::GenericGrayImage, img::AbstractArray{<:Color3}) =
    f(out, of_eltype(Gray, img))
