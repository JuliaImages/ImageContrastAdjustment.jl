
"""
```
    MidwayEqualization <: AbstractHistogramAdjustmentAlgorithm
    MidwayEqualization(; nbins = 256, minval = 0, maxval = 1)

    adjust_histogram([T,] img_sequence, f::MidwayEqualization(nbins = 256, edges = nothing))
    adjust_histogram!([out_sequence,] img_sequence, f::MidwayEqualization(nbins = 256, edges = nothing))
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

Various options for the parameters of the `adjust_histogram` function and
`MidwayEqualization` types are described in more detail below.

## Choices for `img_sequence`

The function `adjust_histogram` expects a length-2 `Vector` of images (the pair
of images) and returns a length-2 `Vector` of modified images.  The  function
can handle a variety of input types. The type of the returned image matches the
input type.

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
img1 = adjust_histogram(img, GammaCorrection(gamma = 2))
img2 = adjust_histogram(img, GammaCorrection(gamma = 1.2))

# Midway histogram equalization will transform these two images so that their
# intensity distributions are almost identical.
img_sequence = adjust_histogram([img1, img2], MidwayEqualization(nbins = 256))
img1o = first(img_sequence)
img2o = last(img_sequence)
```

# References
1. T. Guillemot and J. Delon, “*Implementation of the Midway Image Equalization*,” Image Processing On Line, vol. 5, pp. 114–129, Jun. 2016. [doi:10.5201/ipol.2016.140](https://doi.org/10.5201/ipol.2016.140)
"""
Base.@kwdef struct MidwayEqualization{T₁ <: Union{Integer, Nothing},
                					  T₂ <: Union{AbstractRange, Nothing}} <: AbstractHistogramAdjustmentAlgorithm
    nbins::T₁ = 256
    edges::T₂ = nothing
end

function (f::MidwayEqualization)(out_sequence::Vector{<:GenericGrayImage}, in_sequence::Vector{<:GenericGrayImage})
    length(out_sequence) == 2 || error("Please supply a length-2 output vector to store the pair of images.")
    length(in_sequence) == 2 || error("Please supply a length-2 input vector storing the image pair.")
    out1 = first(out_sequence)
    out2 = last(out_sequence)
    in1 = first(in_sequence)
    in2 = last(in_sequence)
    out1 .= in1
    out2 .= in2
    edges, pdf1, pdf2 = isnothing(f.edges) ? construct_pdfs(out1, out2, f.nbins) : construct_pdfs(out1, out2, f.edges)
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
    index₁ = firstindex(out_sequence)
    index₂ = lastindex(out_sequence)
    out_sequence[index₁]  = match_pdf!(out1, edges, pdf1, midway_pdf)
    out_sequence[index₂]  = match_pdf!(out2, edges, pdf2, midway_pdf)
end

function (f::MidwayEqualization)(out_sequence::Vector{<:AbstractArray{<:Color3}}, in_sequence::Vector{<:AbstractArray{<:Color3}})
    length(out_sequence) == 2 || error("Please supply a length-2 output vector to store the pair of images.")
    length(in_sequence) == 2 || error("Please supply a length-2 input vector storing the image pair.")
    out1 = first(out_sequence)
    out2 = last(out_sequence)
    in1 = first(in_sequence)
    in2 = last(in_sequence)

    T₁ = eltype(out1)
    T₂ = eltype(out2)

    out1 .= in1
    out2 .= in2

    in_yiq1 = convert.(YIQ, in1)
    in_yiq1_view = channelview(in_yiq1)
    in_yiq2 = convert.(YIQ, in2)
    in_yiq2_view = channelview(in_yiq2)
    in_yiq_view_sequence = [view(in_yiq1_view,1,:,:), view(in_yiq2_view,1,:,:)]

    out_yiq1 = convert.(YIQ, out1)
    out_yiq1_view = channelview(out_yiq1)
    out_yiq2 = convert.(YIQ, out2)
    out_yiq2_view = channelview(out_yiq2)
    out_yiq_view_sequence = [view(out_yiq1_view,1,:,:), view(out_yiq2_view,1,:,:)]

    adjust_histogram!(out_yiq_view_sequence, in_yiq_view_sequence, f)

    index₁ = firstindex(out_sequence)
    index₂ = lastindex(out_sequence)
    out_sequence[index₁] .= convert.(T₁, out_yiq1)
    out_sequence[index₂] .= convert.(T₂, out_yiq2)
end

(f::MidwayEqualization)(out_sequence::Vector{<:GenericGrayImage}, img_sequence::Vector{<:AbstractArray{<:Color3}}) =
    f(out_sequence, map(img -> of_eltype(Gray, img), img_sequence))
