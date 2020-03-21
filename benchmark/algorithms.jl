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

# # print_markdown locally
# function print_markdown(io, judgement)
#     buffer = IOBuffer()
#     for alg_name in keys(judgement)
#         println(buffer, "# ", alg_name, "\n")
#         println(buffer, "Parameters")
#         println(buffer, "* sizes: ", join(zip(string.(tst_sizes), string.(tst_sizes)), ", "))
#         println(buffer, "* types: ", join(tst_types, ", "))
#         println(buffer)

#         println(buffer, "| size  |  type  |  time  |  memory |")
#         println(buffer, "| -------| ----------| -------| -------- |")
#         for T in keys(judgement[alg_name])
#             for sz in keys(judgement[alg_name][T])
#                 # t = judgement["LinearStretching_custom"][Float32]["64×64"]
#                 t = judgement[alg_name][T][sz]
#                 _time = @sprintf "%.6f" t.time * 1e-6 # ms
#                 _memory = @sprintf "%.6f" t.memory * 1e-6 # MiB
#                 println(buffer, "| ", T,
#                             " | ", sz,
#                             " | ", _time, " ms",
#                             " | ", _memory, " MiB",
#                             "|")
#             end
#         end
#         println(buffer, "\n")
#     end
#     print(io, String(take!(buffer)))
# end
