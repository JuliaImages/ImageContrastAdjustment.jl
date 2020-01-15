using ImageContrastAdjustment
using Test, ImageCore, ImageFiltering, TestImages, LinearAlgebra

@testset "ImageContrastAdjustment.jl" begin
    include("core.jl")
    include("adaptive_equalization.jl")
    include("histogram_construction.jl")
    include("histogram_matching.jl")
    include("histogram_equalization.jl")
    include("histogram_midway_equalization.jl")
    include("gamma_adjustment.jl")
    include("linear_stretching.jl")
    include("contrast_stretching.jl")
end
