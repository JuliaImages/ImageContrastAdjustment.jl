push!(LOAD_PATH,"../src/")
using Documenter, ImageContrastAdjustment, ColorVectorSpace, ColorTypes
makedocs(sitename="Documentation",
            format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"))
deploydocs(repo = "github.com/JuliaImages/ImageContrastAdjustment.jl.git")
