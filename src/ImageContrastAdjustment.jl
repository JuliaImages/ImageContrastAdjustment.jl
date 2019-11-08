module ImageContrastAdjustment

using ImageCore, ColorTypes, ColorVectorSpace, FixedPointNumbers
using MappedArrays: of_eltype

# TODO: port HistogramAdjustmentAPI to ImagesAPI
include("HistogramAdjustmentAPI/HistogramAdjustmentAPI.jl")
import .HistogramAdjustmentAPI: AbstractHistogramAdjustmentAlgorithm,
                                adjust_histogram, adjust_histogram!

# TODO Relax this to all image color types
const GenericGrayImage = AbstractArray{<:Union{Number, AbstractGray}}

abstract type AbstractHistogramOperation <: AbstractHistogramAdjustmentAlgorithm  end
# TODO Depracate AbstractHistogramOperation and replace with AbstractHistogramAdjustmentAlgorithm
Base.@kwdef struct Equalization{T₁ <: Union{Real,AbstractGray},
								T₂ <: Union{Real,AbstractGray}} <: AbstractHistogramOperation
	nbins::Int = 256
	minval::T₁ = 0.0
	maxval::T₂ = 1.0
end

Base.@kwdef struct MidwayEqualization{T₁ <: Union{Integer, Nothing},
                					   T₂ <: Union{AbstractRange, Nothing}} <: AbstractHistogramOperation
	nbins::T₁ = 256
	edges::T₂ = nothing
end

Base.@kwdef struct Matching{T₁ <: AbstractArray,
                T₂ <: Union{Integer, Nothing},
                T₃ <: Union{AbstractRange, Nothing}} <: AbstractHistogramOperation
    targetimg::T₁
    nbins::T₂ = 256
    edges::T₃ = nothing
end

Base.@kwdef struct GammaCorrection{T <: Real} <: AbstractHistogramOperation
	gamma::T = 1.0
end

Base.@kwdef struct LinearStretching{T₁ <: Union{Real,AbstractGray},
							        T₂ <: Union{Real,AbstractGray}} <: AbstractHistogramOperation
	minval::T₁ = 0.0
	maxval::T₂ = 1.0
end

Base.@kwdef struct ContrastStretching{T₁ <: Union{Real,AbstractGray},
						              T₂ <: Union{Real,AbstractGray}}  <: AbstractHistogramOperation
     t::T₁ = 0.5
     slope::T₂ = 1.0
end


include("core.jl")
include("build_histogram.jl")
include("algorithms/common.jl")
include("algorithms/equalization.jl")
include("algorithms/linear_stretching.jl")
include("algorithms/contrast_stretching.jl")
include("algorithms/gamma_correction.jl")
include("algorithms/matching.jl")
include("algorithms/midway_equalization.jl")
#include("contrastadjustment.jl")

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
