using BenchmarkTools
using ImageContrastAdjustment
using ImageTransformations
using TestImages
using ImageCore
using Logging, TerminalLoggers
using Printf

include("algorithms.jl")

on_CI = haskey(ENV, "GITHUB_ACTIONS")

const SUITE = BenchmarkGroup()

img = testimage("mandril_gray")
alg_list = (( "LinearStretching", LinearStretching()),
            ( "AdaptiveEqualization", AdaptiveEqualization()),
            ( "ContrastStretching", ContrastStretching()),
            ( "Equalization", Equalization()),
            ( "GammaCorrection", GammaCorrection()),
            ( "Matching", Matching),
            ( "MidwayEqualization", MidwayEqualization()))

tst_sizes = on_CI ? (64, 256, ) : (64, 128, 256, 512)
tst_types = (N0f8, Float32, Gray{N0f8}, Gray{Float32}, RGB{N0f8}, RGB{Float32})

for (alg_name, alg) in alg_list
    add_algorithm_benchmark!(SUITE, img, alg_name, alg;
                             tst_sizes=tst_sizes,
                             tst_types=tst_types)
end


# # run benchmark locally
# if !on_CI
#     with_logger(TerminalLogger()) do

#         if !isfile(@__DIR__, "params.json")
#             tuning = tune!(SUITE; verbose = true);
#             BenchmarkTools.save("params.json", "SUITE", params(SUITE))
#         end
#         loadparams!(SUITE, BenchmarkTools.load("params.json")[1], :evals, :samples);

#         results = run(SUITE; verbose=true)
#         judgement = median(results)

#         if get(ENV, "CI", "local") == "local"
#             open(joinpath(@__DIR__, "reports.md"), "w") do io
#                 print_markdown(io, judgement)
#             end
#         end
#     end
# end
