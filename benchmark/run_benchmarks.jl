# To run it locally, BenchmarkCI and JSON should be added to root project
using BenchmarkCI
using JSON
on_CI = haskey(ENV, "GITHUB_ACTIONS")

# github action workflow doesn't support action on labels
# this function manually check if PR has label "run benchmarks"
function should_do_benchmark(event_path = get(ENV, "GITHUB_EVENT_PATH", nothing))
    event_path === nothing && return false
    event = JSON.parsefile(event_path)
    haskey(event, "pull_request") || return false
    labels = [x["name"] for x in event["pull_request"]["labels"]]
    return "run benchmarks" in labels
end

do_benchmark = on_CI ? should_do_benchmark() : true

if do_benchmark
    BenchmarkCI.judge()
    on_CI ? BenchmarkCI.postjudge() : BenchmarkCI.displayjudgement()
end
