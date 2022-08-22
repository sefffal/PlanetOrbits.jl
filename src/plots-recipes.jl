#=
This file contains plot recipes for the Plots.jl
ecosystem. This way you can do e.g.:
plot(elems)
=#

# Plotting recipes for orbital elements
using RecipesBase
@recipe function f(elem::AbstractOrbit)
    os = orbitsolve_ν(elem, 0)
    line_z --> nothing
    solmarker --> false
    os
end

# Recipe for an array of orbits. Same as sigle orbit,
# but scale down transparency as we add more orbits.  
@recipe function f(elems::AbstractArray{<:AbstractOrbit})
    label --> ""
    seriesalpha --> 30/length(elems)
    for elem in elems
        @series begin
            elem
        end
    end
end


default_plotkind(os::OrbitSolutionKep) = (:x, :y)
default_plotkind(os::OrbitSolutionVisual) = :astrometry
default_plotkind(os::OrbitSolutionRadialVelocity) = :radvel

# Plotting recipes for orbital elements
using RecipesBase
@recipe function f(os::AbstractOrbitSolution)
    
    if isdefined(Main, :Plots) && isdefined(Main.Plots, :palette) && get(plotattributes, :solmarker, true) 
        # We would like to create a nice semi-transparent 
        # gray gradient. But palette isn't in PlotRecipes so
        # we fall back to this hacky way of getting it
        color --> Main.Plots.palette(["#444444ff", "#44444433"],10)
    end

    kind = get(plotattributes, :kind, default_plotkind(os))

    # We trace out in equal steps of true anomaly instead of time for a smooth
    # curve, regardless of eccentricity.
    eccanoms = range(os.EA, os.EA+2π, length=90)
    solns = orbitsolve_eccanom.(os.elem, eccanoms)

    line_z --> eccanoms
    colorbar --> nothing

    if kind == :astrometry
        xs = raoff.(solns)
        ys = decoff.(solns)

        # We almost always want to see spatial coordinates with equal step sizes
        aspect_ratio --> 1
        # And we almost always want to reverse the RA coordinate to match how we
        # see it in the sky.
        xflip --> true

        xguide --> "Δra - (mas)"
        yguide --> "Δdec - (mas)"

        @series begin
            xs, ys
        end
        if get(plotattributes, :solmarker, true)
            @series begin
                seriestype --> :scatter
                label --> ""
                color --> :gray
                [raoff(os)], [decoff(os)]
            end
        end

    elseif kind isa Tuple || kind ∈ (:x, :y, :z)
        if kind isa Symbol
            kind = (kind,)
        end
        varx = kind[1]
        if length(kind) > 1
            vary = kind[2]
        end
        if length(kind) > 2
            varz = kind[3]
        end
        funcs = (;
            x=PlanetOrbits.posx,
            y=PlanetOrbits.posy,
            z=PlanetOrbits.posz,
        )
        guides = (;
            x="x (au)",
            y="y (au)",
            z="z (au)",
        )

        xs = funcs[varx].(solns)
        x = xs[begin]
        if length(kind) > 1
            ys = funcs[vary].(solns)
            y = ys[begin]
        end
        if length(kind) > 2
            zs = funcs[varz].(solns)
            z = zs[begin]
        end


        # Set the guides if unset
        if length(kind) == 1
            xguide --> "time (days)"
            yguide --> guides[varx]
        else
            xguide --> guides[varx]
            # We almost always want to see spatial coordinates with equal step sizes
            aspect_ratio --> 1
            xflip --> true
        end
        if length(kind) > 1
            yguide --> guides[vary]
        end
        if length(kind) > 2
            zguide --> guides[varz]
        end

        if length(kind) == 1
            ts = _time_from_soln.(solns)
            t = ts[begin]
            # Prevent wrapped line
            i = findfirst(ts[begin:end-1] .> ts[begin+1:end])
            xs[i] = NaN
            @series begin
                ts, xs
            end
            if get(plotattributes, :solmarker, true)
                @series begin
                    seriestype --> :scatter
                    label --> ""
                    [t], [x]
                end
            end
        elseif length(kind) == 2
            @series begin
                xs, ys
            end
            if get(plotattributes, :solmarker, true)
                @series begin
                    seriestype --> :scatter
                    label --> ""
                    color --> :gray
                    [x], [y]
                end
            end
        elseif length(kind) == 3
            @series begin
                xs, ys, zs
            end
            if get(plotattributes, :solmarker, true)
                @series begin
                    seriestype --> :scatter
                    label --> ""
                    color --> :gray
                    [x],[y],[z]
                end
            end
        else
            error("Cannot plot in more than 3 dimensions")
        end
    elseif kind == :radvel
        rvs = radvel.(solns)
        ts = _time_from_soln.(solns)
        # Prevent wrapped line
        i = findfirst(ts[begin:end-1] .> ts[begin+1:end])
        rvs[i] = NaN

        xguide --> "time (days)"
        yguide --> "secondary radial velocity (m/s)"

        @series begin
            ts, rvs
        end
        if get(plotattributes, :solmarker, true)
            @series begin
                seriestype --> :scatter
                label --> ""
                color --> :gray
                ts[begin:begin], rvs[begin:begin]
            end
        end
    else
        error("unrecognized kind=$kind value")
    end


end
@recipe function f(oses::AbstractArray{<:OrbitSolutionVisual})

    label --> ""
    seriesalpha --> 30/length(oses)
    for os in oses
        @series begin
            os
        end
    end
end


function _time_from_soln(os::AbstractOrbitSolution; tref=58849)

    # ---- Worked math --- 
    # ν/2 = atan(elem.ν_fact*tan(EA/2))
    # tan(ν/2) = elem.ν_fact*tan(EA/2)
    # tan(ν/2)/elem.ν_fact = tan(EA/2)
    # atan(tan(ν/2)/elem.ν_fact) = (EA/2)
    # atan(tan(ν/2)/elem.ν_fact)*2 = EA
    # EA = atan(tan(ν/2)/elem.ν_fact)*2

    # Compute eccentric anomaly
    MA = os.EA - os.elem.e * sin(os.EA)

    # Epoch of periastron passage
    tₚ = periastron(os.elem, tref)
    
    
    # MA = meanmotion(elem)/oftype(t, year2day) * (t - tₚ)
    # MA /  meanmotion(elem) * year2day = (t - tₚ)
    # MA /  meanmotion(elem) * year2day + tₚ = t
    t = MA /  meanmotion(os.elem) * year2day + tₚ - tref

    return t
end