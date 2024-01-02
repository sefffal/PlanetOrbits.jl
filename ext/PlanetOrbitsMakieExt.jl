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
    νs = range(-π, π, length=90)
    return map(kep2cart_ν.(elem, νs)) do coord
        return Makie.Point2f(coord.x, coord.y)
    end
end

function Makie.convert_single_argument(elems::Vector{<:AbstractOrbit})
    # We trace out in equal steps of true anomaly instead of time for a smooth
    # curve, regardless of eccentricity.
    νs = range(-π, π, length=90)

    coords = kep2cart_ν.(elems, νs')

    xs = [c[1] for c in coords]'
    ys = [c[2] for c in coords]'

    # Treat as one long series interrupted by NaN
    xs = reduce(vcat, [[x; NaN] for x in eachcol(xs)])
    ys = reduce(vcat, [[y; NaN] for y in eachcol(ys)])

    return map(zip(xs,ys)) do (x,y)
        return Makie.Point2f(x, y)
    end
end

end