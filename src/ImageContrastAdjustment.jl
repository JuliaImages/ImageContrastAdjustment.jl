module ImageContrastAdjustment

using ImageCore, ColorTypes, ColorVectorSpace, FixedPointNumbers
using MappedArrays: of_eltype

# TODO: port HistogramAdjustmentAPI to ImagesAPI
include("HistogramAdjustmentAPI/HistogramAdjustmentAPI.jl")
import .HistogramAdjustmentAPI: AbstractHistogramAdjustmentAlgorithm,
                                adjust_histogram, adjust_histogram!

# TODO Relax this to all image color types
const GenericGrayImage = AbstractArray{<:Union{Number, AbstractGray}}

include("core.jl")
include("build_histogram.jl")
include("algorithms/common.jl")
include("algorithms/equalization.jl")
include("algorithms/linear_stretching.jl")
include("algorithms/contrast_stretching.jl")
include("algorithms/gamma_correction.jl")
include("algorithms/matching.jl")
include("algorithms/midway_equalization.jl")

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
