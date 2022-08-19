using Documenter, PlanetOrbits

include("pages.jl")
makedocs(
    sitename="PlanetOrbits.jl",
    pages=pages,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)


deploydocs(
    repo = "github.com/sefffal/PlanetOrbits.jl.git",
    devbranch = "master"
)
