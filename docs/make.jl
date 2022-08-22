using Documenter, PlanetOrbits

ENV["GKSwstype"] = "100" # GR in documeter env variable

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
