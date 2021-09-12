module ImageContrastAdjustment

using ImageCore
using ImageTransformations: imresize
# Where possible we avoid a direct dependency to reduce the number of [compat] bounds
using ImageCore.MappedArrays
using Parameters: @with_kw # Same as Base.@kwdef but works on Julia 1.0

# TODO: port HistogramAdjustmentAPI to ImagesAPI
include("ContrastAdjustmentAPI/ContrastAdjustmentAPI.jl")
import .ContrastAdjustmentAPI: AbstractHistogramAdjustmentAlgorithm,
                                adjust_histogram, adjust_histogram!,
                              AbstractIntensityAdjustmentAlgorithm,
                                adjust_intensity, adjust_intensity!

# TODO Relax this to all image color types
const GenericGrayImage = AbstractArray{<:Union{Number, AbstractGray}}

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
include("compat.jl")
include("deprecations.jl")

export
    # main types and functions
    AdaptiveEqualization,
    Equalization,
    MidwayEqualization,
    Matching,
    GammaCorrection,
    LinearStretching,
    ContrastStretching,
    build_histogram,
    adjust_histogram,
    adjust_histogram!,
    adjust_intensity,
    adjust_intensity!

end # module
