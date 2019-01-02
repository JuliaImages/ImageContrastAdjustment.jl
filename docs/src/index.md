# ImageContrastAdjustment.jl Documentation

## Functions

```@docs
build_histogram
adjust_histogram(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})
adjust_histogram!(::Equalization, ::AbstractArray, ::Integer; ::Union{Real,AbstractGray}, ::Union{Real,AbstractGray})
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
