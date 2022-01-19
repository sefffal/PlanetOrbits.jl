#=
This file contains plot recipes for the Makie.jl
ecosystem. This way you can do e.g.:
lines(elems)
=#

function Makie.convert_single_argument(elem::KeplerianElements)
    # We trace out in equal steps of true anomaly instead of time for a smooth
    # curve, regardless of eccentricity.
    νs = range(-π, π, length=90)
    return map(kep2cart_ν.(elem, νs)) do coord
        return Makie.Point2f(coord.x, coord.y)
    end
end