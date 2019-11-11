
"""
```
    LinearStretching <: AbstractHistogramAdjustmentAlgorithm
    LinearStretching(; minval = 0, maxval = 1)

    adjust_histogram([T,] img, f::LinearStretching)
    adjust_histogram!([out,] img, f::LinearStretching)
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

Various options for the parameters of the `adjust_histogram` and
`LinearStretching` type  are described in more detail below.

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
using ImageContrastAdjustment, ImageView, TestImages

img = testimage("mandril_gray")
imgo = adjust_histogram(img, LinearStretching(minval = 0, maxval = 1))

```

# References
1. W. Burger and M. J. Burge. *Digital Image Processing*. Texts in Computer Science, 2016. [doi:10.1007/978-1-4471-6684-9](https://doi.org/10.1007/978-1-4471-6684-9)

"""
Base.@kwdef struct LinearStretching{T₁ <: Union{Real,AbstractGray},
                                    T₂ <: Union{Real,AbstractGray}} <: AbstractHistogramAdjustmentAlgorithm
    minval::T₁ = 0.0
    maxval::T₂ = 1.0
end

function (f::LinearStretching)(out::GenericGrayImage, img::GenericGrayImage)
    src_minval = minfinite(img)
    src_maxval = maxfinite(img)
    T = eltype(out)
    out .= img
    map!(out,out) do val
        if isnan(val)
            return val
        else
            newval = linear_stretch(val, src_minval, src_maxval, f.minval, f.maxval)
            return  T <: Integer ? round(Int, newval ) : newval
        end
    end
end

function (f::LinearStretching)(out::AbstractArray{<:Color3}, img::AbstractArray{<:Color3})
    T = eltype(out)
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    adjust_histogram!(view(yiq_view,1,:,:), f)
    out .= convert.(T, yiq)
end

(f::LinearStretching)(out::GenericGrayImage, img::AbstractArray{<:Color3}) =
    f(out, of_eltype(Gray, img))
