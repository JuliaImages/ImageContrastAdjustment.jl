# ImageContrastAdjustment.jl Documentation

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

## Functions

```@docs
build_histogram
adjust_histogram(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})
adjust_histogram!(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})
adjust_histogram(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer)
adjust_histogram!(::MidwayEqualization, ::AbstractArray, ::AbstractArray, ::Integer)
adjust_histogram(::Matching, ::AbstractArray, ::AbstractArray, ::Integer)
adjust_histogram!(::Matching, ::AbstractArray, ::AbstractArray, ::Integer)
adjust_histogram(::GammaCorrection, ::AbstractArray, ::Real)
adjust_histogram!(::GammaCorrection, ::AbstractArray, ::Real)
adjust_histogram(::LinearStretching, ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})
adjust_histogram!(::LinearStretching, ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})
adjust_histogram(::ContrastStretching, ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})
adjust_histogram!(::ContrastStretching, ::AbstractArray; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})
```
## Index

```@index
```
