module ImageContrastAdjustment

using ImageCore, ColorTypes, ColorVectorSpace, FixedPointNumbers

abstract type AbstractHistogramOperation end
struct Equalization <: AbstractHistogramOperation end
struct Matching <: AbstractHistogramOperation end

include("core.jl")
include("contrastadjustment.jl")

export
	# main types and functions
	Equalization,
	Matching,
    build_histogram,
	adjust_histogram,
	adjust_histogram!

end # module
