
"""
```
    GammaCorrection <: AbstractIntensityAdjustmentAlgorithm
    GammaCorrection(; gamma = 1)

    adjust_intensity([T,] img, f::GammaCorrection)
    adjust_intensity!([out,] img, f::GammaCorrection)
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

Various options for the parameters of the `adjust_intensity` function and the
`Gamma` type are described in more detail below.

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
imgadj = adjust_intensity( img, GammaCorrection(gamma = 1/2))

# Display the original and adjusted image.
imshow(img)
imshow(imgadj)
```

# References
1. W. Burger and M. J. Burge. *Digital Image Processing*. Texts in Computer Science, 2016. [doi:10.1007/978-1-4471-6684-9](https://doi.org/10.1007/978-1-4471-6684-9)
"""
@with_kw struct GammaCorrection{T <: Real} <: AbstractIntensityAdjustmentAlgorithm
    gamma::T = 1.0
end

function (f::GammaCorrection)(out::GenericGrayImage, img::GenericGrayImage)
    out .= img
    correct_gamma!(out, f.gamma)
    return out
end

# TODO Expand this to more generic gray color types.
function correct_gamma!(img::AbstractArray{Gray{T}},  gamma::Real) where T <: FixedPointNumbers.Normed
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

function correct_gamma!(img::GenericGrayImage,  gamma::Real)
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

function (f::GammaCorrection)(out::AbstractArray{<:Color3}, img::AbstractArray{<:Color3})
    T = eltype(out)
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    adjust_histogram!(view(yiq_view,1,:,:), f)
    out .= convert.(T, yiq)
end

(f::GammaCorrection)(out::GenericGrayImage, img::AbstractArray{<:Color3}) =
    f(out, of_eltype(Gray, img))
