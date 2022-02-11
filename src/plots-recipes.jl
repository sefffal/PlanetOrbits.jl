#=
This file contains plot recipes for the Plots.jl
ecosystem. This way you can do e.g.:
plot(elems)
=#

# Plotting recipes for orbital elementst
using RecipesBase
@recipe function f(elem::AbstractElements)
    # We trace out in equal steps of true anomaly instead of time for a smooth
    # curve, regardless of eccentricity.
    νs = range(-π, π, length=90)
    solns = orbitsolve_ν.(elem, νs)
    xs = raoff.(solns)
    ys = decoff.(solns)

    # We almost always want to see spatial coordinates with equal step sizes
    aspect_ratio --> 1
    # And we almost always want to reverse the RA coordinate to match how we
    # see it in the sky.
    xflip --> true

    return xs, ys
end

# Recipe for an array of orbits. Same as sigle orbit,
# but scale down transparency as we add more orbits.  
@recipe function f(elems::AbstractArray{<:AbstractElements})

    # Step through true anomaly instead of time.
    # This produces far nicer looking plots, especially if
    # the orbits in question vary significantly in period
    # or are eccentric
    νs = range(-π, π, length=90)
    solns = orbitsolve_ν.(elems, νs')

    xs = raoff.(solns)'
    ys = decoff.(solns)'

    # Treat as one long series interrupted by NaN
    xs = reduce(vcat, [[x; NaN] for x in eachcol(xs)])
    ys = reduce(vcat, [[y; NaN] for y in eachcol(ys)])

    # We almost always want to see spatial coordinates with equal step sizes
    aspect_ratio --> 1
    # And we almost always want to reverse the RA coordinate to match how we
    # see it in the sky.
    xflip --> true
    xguide --> "ΔRA - mas"
    yguide --> "ΔDEC - mas"

    seriesalpha --> 30/length(elems)

    return xs, ys
end