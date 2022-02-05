using Documenter, DirectOrbits


makedocs(
    sitename="DirectOrbits.jl",
    pages = [
        "Home" => "index.md",
        # "Getting Started" => "getting-started.md",
        "Tutorials" => [
            "Plotting" => "plots.md",
            "Image Warping" => "image-warping.md",
        ],
        "Documentation" => [
            "Conventions" => "conventions.md",
            "Kepler Solver" => "kepler.md",
            "API" => "api.md"
        ]
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)


deploydocs(
    repo = "github.com/sefffal/DirectOrbits.jl.git",
    devbranch = "master"
)
