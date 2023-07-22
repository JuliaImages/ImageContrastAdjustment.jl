using ImageContrastAdjustment
using Test, ImageCore, ImageFiltering, TestImages, LinearAlgebra
using Aqua

if Base.VERSION >= v"1.6"
    @testset "Aqua" begin
        # TODO: fix the ambiguity test
        Aqua.test_all(ImageContrastAdjustment; ambiguities=false)
    end
end

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
