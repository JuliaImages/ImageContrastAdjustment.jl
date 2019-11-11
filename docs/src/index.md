# ImageContrastAdjustment.jl Documentation

A Julia package for enhancing and manipulating image contrast.

```@contents
Depth = 3
```

## Getting started
This package is part of a wider Julia-based image processing
[ecosystem](https://github.com/JuliaImages). If you are starting out, then you
may benefit from [reading](https://juliaimages.org/latest/quickstart/) about
some fundamental conventions that the ecosystem utilizes that are markedly
different from how images are typically represented in OpenCV, MATLAB, ImageJ or
Python.

The usage examples in the `ImageContrastAdjustment.jl` package assume that you have
already installed some key packages. Notably, the examples assume that you are
able to load and display an image. Loading an image is facilitated through the
[FileIO.jl](https://github.com/JuliaIO/FileIO.jl) package, which uses
[QuartzImageIO.jl](https://github.com/JuliaIO/QuartzImageIO.jl) if you are on
`MacOS`, and [ImageMagick.jl](https://github.com/JuliaIO/ImageMagick.jl)
otherwise. Depending on your particular system configuration, you might
encounter problems installing the image loading packages, in which case you can
refer to the [troubleshooting
guide](https://juliaimages.org/latest/troubleshooting/#Installation-troubleshooting-1).

Image display is typically handled by the
[ImageView.jl](https://github.com/JuliaImages/ImageView.jl) package. However,
there are some known issues with this package. For example, on `Windows` the
package has the side-effect of introducing substantial [input
lag](https://github.com/JuliaImages/ImageView.jl/issues/176) when typing in the
Julia REPL. Also, as of writing, some users of `MacOS` are [unable to
use](https://github.com/JuliaImages/ImageView.jl/issues/175) the `ImageView.jl`
package.

As an alternative, one can display an image using the
[Makie.jl](https://github.com/JuliaPlots/Makie.jl) plotting package. There is
also the [ImageShow.jl](https://github.com/JuliaImages/ImageShow.jl) package
which facilitates displaying images in `Jupyter` notebooks via
[IJulia.jl](https://github.com/JuliaLang/IJulia.jl).

Finally, one can also obtain a useful preview of an image in the REPL using the
[ImageInTerminal.jl](https://github.com/JuliaImages/ImageInTerminal.jl) package.
However, this package assumes that the terminal uses a monospace font, and tends
not to produce adequate results in a Windows environment.

Another package that is used to illustrate the functionality in
`ImageContrastAdjustment.jl` is the
[TestImages.jl](https://github.com/JuliaImages/TestImages.jl) which serves as a
repository of many standard image processing test images.


## Basic usage

Each contrast manipulation algorithm in `ImageContrastAdjustment.jl` is an
[`AbstractHistogramAdjustmentAlgorithm`](@ref
ImageContrastAdjustment.HistogramAdjustmentAPI.AbstractHistogramAdjustmentAlgorithm).

Suppose one wants to enhance the contrast an image. This can be achieved by
simply choosing an appropriate algorithm and calling [`adjust_histogram`](@ref)
or [`adjust_histogram!`](@ref) in the image. The contrast will be automatically
enhanced.

Let's see a simple demo:

```@example
using TestImages, ImageContrastAdjustment
using FileIO # hide
img = testimage("cameraman")
alg = Equalization(nbins = 256)
img_adjusted = adjust_histogram(img, alg)
save("images/demo.jpg", hcat(img, img_adjusted)) # hide
```

```@raw html
<img src="images/demo.jpg" width="400px" alt="demo image" />
```

This usage reads as "`adjust_histogram` of the image `img` with algorithm `alg`"

For more advanced usage, please check [function reference](@ref function_reference) page.

## Examples
Below are some examples of contrast adjustment operations that this package facilitates.

```@raw html
<table width="500" border="0" cellpadding="5">

<tr>

<td align="center" valign="center">
<img src="images/contrast_stretching.gif" width="100px" alt="Contrast Stretching" />
<br />
Contrast Stretching
</td>

<td align="center" valign="center">
<img src="images/contrast_stretching_col.gif" width="100px" alt="Contrast Stretching" />
<br />
Contrast Stretching
</td>

</tr>

<tr>

<td align="center" valign="center">
<img src="images/linear_stretching.gif" width="100px" alt="Linear Stretching (Normalization)" />
<br />
Linear Stretching/Normalization
</td>

<td align="center" valign="center">
<img src="images/linear_stretching_col.gif" width="100px" alt="Linear Stretching (Normalization)" />
<br />
Linear Stretching/Normalization
</td>

</tr>

<tr>

<td align="center" valign="center">
<img src="images/gamma_correction.gif" width="100px" alt="Gamma Correction" />
<br />
Gamma Correction
</td>

<td align="center" valign="center">
<img src="images/gamma_correction_col.gif" width="100px" alt="Gamma Correction" />
<br />
Gamma Correction
</td>

</tr>

<tr>

<td align="center" valign="center">
<img src="images/equalization.gif" width="100px" alt="Histogram Equalization" />
<br />
Histogram Equalization
</td>

<td align="center" valign="center">
<img src="images/equalization_col.gif" width="100px" alt="Histogram Equalization" />
<br />
Histogram Equalization
</td>

</tr>

<tr>

<td align="center" valign="center">
<img src="images/midway_equalization.gif" width="100px" alt="Midway Histogram Equalization" />
<br />
Midway Histogram Equalization
</td>

<td align="center" valign="center">
<img src="images/midway_equalization_col.gif" width="100px" alt="Midway Histogram Equalization" />
<br />
Midway Histogram Equalization
</td>

</tr>

<tr>

<td align="center" valign="center">
<img src="images/matching.gif" width="100px" alt="Histogram Matching" />
<br />
Histogram Matching
</td>

<td align="center" valign="center">
<img src="images/matching_col.gif" width="100px" alt="Histogram Matching" />
<br />
Histogram Matching
</td>

</tr>

</table>
```
