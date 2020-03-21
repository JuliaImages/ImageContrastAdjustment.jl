# Usage:
#     julia benchmark/run_benchmarks.jl
using BenchmarkTools
using ImageContrastAdjustment
using ImageTransformations
using TestImages
using ImageCore

on_CI = haskey(ENV, "GITHUB_ACTIONS")

img = testimage("mandril_gray")
tst_sizes = on_CI ? (64, ) : (64, 128)
tst_types = (N0f8, Float32, Gray{N0f8}, Gray{Float32}, RGB{N0f8}, RGB{Float32})

const SUITE = BenchmarkGroup()

# TODO: add MidwayEqualization
alg_list = (( "LinearStretching", LinearStretching()),
            ( "AdaptiveEqualization", AdaptiveEqualization()),
            ( "ContrastStretching", ContrastStretching()),
            ( "Equalization", Equalization()),
            ( "GammaCorrection", GammaCorrection()),
            ( "Matching", Matching))


function add_algorithm_benchmark!(suite, img, alg_name, alg;
                                  tst_sizes, tst_types)
    haskey(suite, alg_name) || (suite[alg_name] = BenchmarkGroup())

    for T in tst_types
        haskey(suite[alg_name], T) || (suite[alg_name][T] = BenchmarkGroup())
        for sz in tst_sizes
            tst_img = imresize(T.(img), (sz, sz))
            if alg === Matching
                tst_alg = alg(targetimg = tst_img)
            else
                tst_alg = alg
            end
            
            suite[alg_name][T]["$sz×$sz"] = @benchmarkable adjust_histogram($tst_img, $tst_alg)
        end
    end
end


for (alg_name, alg) in alg_list
    add_algorithm_benchmark!(SUITE, img, alg_name, alg;
                             tst_sizes=tst_sizes,
                             tst_types=tst_types)
end
