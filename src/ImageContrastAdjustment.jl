module ImageContrastAdjustment

using ImageCore, ColorTypes, ColorVectorSpace, FixedPointNumbers

abstract type AbstractHistogramOperation end
struct Equalization <: AbstractHistogramOperation end
struct MidwayEqualization <: AbstractHistogramOperation end
struct Matching <: AbstractHistogramOperation end
struct GammaCorrection <: AbstractHistogramOperation end
struct LinearStretching <: AbstractHistogramOperation end
struct ContrastStretching <: AbstractHistogramOperation end


include("core.jl")
include("contrastadjustment.jl")

export
	# main types and functions
	Equalization,
	MidwayEqualization,
	Matching,
	GammaCorrection,
	LinearStretching,
	ContrastStretching,
    build_histogram,
	adjust_histogram,
	adjust_histogram!

end # module
