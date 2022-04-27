#=
This file contains plot recipes for the Plots.jl
ecosystem. This way you can do e.g.:
plot(elems)
=#

# Plotting recipes for orbital elements
using RecipesBase
@recipe function f(elem::KeplerianElements)

    kind = get(plotattributes, :kind, :astrometry)

    if kind == :astrometry
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
    elseif kind == :radvel
        # We trace out in equal steps of true anomaly instead of time for a smooth
        # curve, regardless of eccentricity.
        νs = range(-π, π, length=90)
        solns = orbitsolve_ν.(elem, νs)
        rvs = radvel.(solns)
        ts = _time_from_trueanom.(elem, νs)

        xguide --> "time (days)"
        yguide --> "secondary radial velocity (m/s)"

        return ts, rvs
    end
end

# Recipe for an array of orbits. Same as sigle orbit,
# but scale down transparency as we add more orbits.  
@recipe function f(elems::AbstractArray{<:KeplerianElements})

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
    label --> ""

    seriesalpha --> 30/length(elems)

    return xs, ys
end


# Plotting recipes for orbital elementst
using RecipesBase
@recipe function f(os::OrbitSolutionKeplerian)
    
    @series begin
        # Hacky
        if isdefined(Main, :Plots)
            color := Main.Plots.palette(["#444444ff", "#44444433"],10)
        end
        
        # We trace out in equal steps of true anomaly instead of time for a smooth
        # curve, regardless of eccentricity.
        νs = range(os.ν, os.ν+2π, length=90)
        solns = orbitsolve_ν.(os.elem, νs)
        xs = raoff.(solns)
        ys = decoff.(solns)
        line_z := νs
        xs,ys
    end

    # We almost always want to see spatial coordinates with equal step sizes
    aspect_ratio --> 1
    # And we almost always want to reverse the RA coordinate to match how we
    # see it in the sky.
    xflip --> true
    colorbar --> nothing


    @series begin
        color --> :gray
        seriestype --> :scatter
        label --> ""

        [raoff(os)], [decoff(os)]
    end
end
@recipe function f(oses::AbstractArray{<:OrbitSolutionKeplerian})

    label --> ""
    seriesalpha --> 30/length(oses)
    for os in oses
        @series begin
            os
        end
    end
end




# Plotting recipes for orbital elements
@recipe function f(elem::RadialVelocityElements)
    # We trace out in equal steps of true anomaly instead of time for a smooth
    # curve, regardless of eccentricity.
    νs = range(-π, π, length=90)
    solns = orbitsolve_ν.(elem, νs)
    rvs = radvel.(solns)
    ts = _time_from_trueanom.(elem, νs)

    xguide --> "time (days)"
    yguide --> "secondary radial velocity (m/s)"

    return ts, rvs
end

function _time_from_trueanom(elem::AbstractOrbit, ν; tref=58849)

    # ----
    # ν/2 = atan(elem.ν_fact*tan(EA/2))
    # tan(ν/2) = elem.ν_fact*tan(EA/2)
    # tan(ν/2)/elem.ν_fact = tan(EA/2)
    # atan(tan(ν/2)/elem.ν_fact) = (EA/2)
    # atan(tan(ν/2)/elem.ν_fact)*2 = EA
    EA = atan(tan(ν/2)/elem.ν_fact)*2

    # Compute eccentric anomaly
    # EA = kepler_solver(MA, elem.e)
    MA = EA - elem.e * sin(EA)

    # Epoch of periastron passage
    tₚ = periastron(elem, tref)
    
    
    # MA = meanmotion(elem)/oftype(t, year2day) * (t - tₚ)
    # MA /  meanmotion(elem) * year2day = (t - tₚ)
    # MA /  meanmotion(elem) * year2day + tₚ = t
    t = MA /  meanmotion(elem) * year2day + tₚ

    return t
end