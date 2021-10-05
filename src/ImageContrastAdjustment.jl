module ImageContrastAdjustment

using ImageCore
using ImageTransformations: imresize
# Where possible we avoid a direct dependency to reduce the number of [compat] bounds
using ImageCore.MappedArrays
using IntervalSets
using Parameters: @with_kw # Same as Base.@kwdef but works on Julia 1.0
using Reexport
using StatsBase: percentile
@reexport using IntervalSets

# TODO: port HistogramAdjustmentAPI to ImagesAPI
include("HistogramAdjustmentAPI/HistogramAdjustmentAPI.jl")
import .HistogramAdjustmentAPI: AbstractHistogramAdjustmentAlgorithm,
                                adjust_histogram, adjust_histogram!

# TODO Relax this to all image color types
const GenericGrayImage = AbstractArray{<:Union{Number, AbstractGray}}

# Temporary definition which will be moved to ImageCore.
"""
    Percentiles(p₁, p₂, p₃, ..., pₙ)
Specifies a length-n list of [percentiles](https://en.wikipedia.org/wiki/Percentile).
"""
struct Percentiles{T} <: Real
    p::T
end
Percentiles(args...) = Percentiles(tuple(args...))
"""
    MinMax()
A type that can be used to signal to [`PiecewiseLinearStretching`](@ref) that it should
use the mininum and maximum value in an image as it source or destination interval. 
"""
struct MinMax end


include("core.jl")
include("build_histogram.jl")
include("algorithms/common.jl")
include("algorithms/adaptive_equalization.jl")
include("algorithms/equalization.jl")
include("algorithms/linear_stretching.jl")
include("algorithms/contrast_stretching.jl")
include("algorithms/gamma_correction.jl")
include("algorithms/matching.jl")
include("algorithms/midway_equalization.jl")
include("algorithms/piecewise_linear_stretching.jl")
include("compat.jl")

export
    # main types and functions
    AdaptiveEqualization,
    Equalization,
    MidwayEqualization,
    Matching,
    GammaCorrection,
    LinearStretching,
    ContrastStretching,
    PiecewiseLinearStretching,
    Percentiles,
    MinMax,
    build_histogram,
    adjust_histogram,
    adjust_histogram!

end # module
