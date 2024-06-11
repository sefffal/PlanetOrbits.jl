#=
This file contains plot recipes for the Makie.jl
ecosystem. This way you can do e.g.:
lines(elems)
=#

module PlanetOrbitsMakieExt

using PlanetOrbits, Makie

function Makie.convert_single_argument(elem::AbstractOrbit)
    # We trace out in equal steps of true anomaly instead of time for a smooth
    # curve, regardless of eccentricity.
    L = 90
    eccanoms = range(-2π, 2π, length=L)
    solns = orbitsolve_eccanom.(elem, eccanoms)
    νs = range(-π, π, length=90)
    return map(solns) do sol
        return Makie.Point2f(raoff(sol), decoff(sol))
    end
end

end