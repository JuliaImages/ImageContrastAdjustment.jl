push!(LOAD_PATH,"../src/")
using Documenter, ImageContrastAdjustment, ColorVectorSpace, ColorTypes
makedocs(sitename="Documentation",
            Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"))
deploydocs(repo = "github.com/zygmuntszpak/ImageContrastAdjustment.jl.git")
