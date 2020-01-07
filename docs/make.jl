using Documenter, ImageContrastAdjustment
makedocs(sitename="Documentation",
            format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"))
deploydocs(repo = "github.com/JuliaImages/ImageContrastAdjustment.jl.git")
