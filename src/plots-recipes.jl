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

    kind = get(plotattributes, :kind, default_plotkind(os))

    if kind == :astrometry
        kind = (:raoff, :decoff)
    end
    if kind isa Symbol
        kind = (:t, kind)
    end

    if length(kind) < 2
        error("Requires at least two variables to create plot")
    end
    

    posvars = (:x,:y,:z)
    astromvars = (:raoff, :decoff)
    astromvelvars = (:pmra, :pmdec)
    ravars = (:x, :raoff, :pmra)
    timevars = (:t, :ν, :trueanom, :meananom, :eccanom)
    spatial_vars = (posvars..., astromvars..., astromvelvars...)
    if all(∈(spatial_vars), kind)
        # Set equal aspect ratio when plotting all spatial coordinates
        aspect_ratio --> 1
    end

    resolver = (;
        t=("t", "days", _time_from_soln,),
        ν=("ν", "rad", trueanom,),
        trueanom=("ν", "rad", trueanom,),
        meananom=("mean.anom.", "rad", meananom,),
        eccanom=("ecc.anom", "rad", eccanom,),
        x=("x", "au", posx,),
        y=("y", "au", posy,), 
        z=("z", "au", posz,),
        raoff=("Δra", "mas", raoff,),
        decoff=("Δdec", "mas", decoff,),
        pmra=("∂ra/∂t", "mas/yr", pmra,),
        pmdec=("∂ra/∂t", "mas/yr", pmdec,),
        accra=("∂²ra/∂t²", "mas/yr²", pmra,),
        accdec=("∂²ra/∂t²", "mas/yr²", pmdec,),
        radvel=("radvel", "m/s", radvel,),
        posangle=("posangle", "mas", posangle,),
        projectedseparation=("projectedseparation", "mas", projectedseparation,),
    )

    xl, xu, xf = resolver[kind[1]]
    yl, yu, yf = resolver[kind[2]]
    xguide --> "$xl ($xu)"
    yguide --> "$yl ($yu)"
    if length(kind) >= 3
        zl, zu, zf = resolver[kind[3]]
        zguide --> "$zl ($zu)"
    end

    bodies = get(plotattributes, :body, :secondary)
    if bodies isa Symbol
        bodies = (bodies,)
    end


    for body in bodies

        
        # We trace out in equal steps of true anomaly instead of time for a smooth
        # curve, regardless of eccentricity.
        # When the independent variable is a timevar (angle or time) we want
        # two cycles, otherwise we just do one complete orbit
        if kind[1] ∈ timevars
            L = 90
            eccanoms = sort([range(-2π, 2π, length=L); os.EA])
            isol = findfirst(==(os.EA), eccanoms)
            solns = orbitsolve_eccanom.(os.elem, eccanoms)
        else
            L = 90
            isol = 1
            eccanoms = range(os.EA, os.EA+2π, length=L)
            solns = orbitsolve_eccanom.(os.elem, eccanoms)
            line_z --> eccanoms
            colorbar --> nothing
        end


        if body == :secondary
            xs = xf.(solns)
            ys = yf.(solns)
            if length(kind) >= 3
                zs = zf.(solns)
            end
        elseif body == :primary
            xs = xf.(solns, plotattributes[:mass])
            ys = yf.(solns, plotattributes[:mass])
            if length(kind) >= 3
                zs = zf.(solns, plotattributes[:mass])
            end
        else
            error("Unrecognized body $body. Pass body=:primary or :secondary")
        end

        # We almost always want to reverse the RA coordinate to match how we
        # see it in the sky.
        if kind[1] ∈ ravars
            xflip --> true
        end
        if kind[2] ∈ ravars
            yflip --> true
        end
        if length(kind) >= 3 && kind[3] ∈ ravars
            zflip --> true
        end

        # Prevent wrapped lines
        if kind[1] ∈ timevars
            P = kind[1]==:t ? period(os.elem) : 2π
            unwrap!(xs, P)
            xs .-= P
        end
        if kind[2] ∈ timevars
            P = kind[2]==:t ? period(os.elem) : 2π
            unwrap!(ys, P)
            ys .-= P
        end
        if length(kind) >= 3 && kind[3] ∈ timevars
            P = kind[3]==:t ? period(os.elem) : 2π
            unwrap!(zs, P)
            zs .-= P
        end

        @series begin
            label --> string(body)
            if isdefined(Main, :Plots) && isdefined(Main.Plots, :palette) && get(plotattributes, :solmarker, true) 
                # We would like to create a nice semi-transparent 
                # gray gradient. But palette isn't in PlotRecipes so
                # we fall back to this hacky way of getting it
                if body == :secondary
                    color --> Main.Plots.palette(["#444444ff", "#44444433"],10)
                else
                    color --> Main.Plots.palette(["#BB4444ff", "#BB444433"],10)
                end
            end
            if length(kind) >= 3
                xs, ys, zs
            else
                xs, ys
            end
        end
        if get(plotattributes, :solmarker, true)
            @series begin
                seriestype --> :scatter
                label --> ""
                if body == :secondary
                    color --> :gray
                else
                    color --> "#BB4444"
                end

                if length(kind) >= 3
                    xs[isol:isol], ys[isol:isol], zs[isol:isol]
                else
                    xs[isol:isol], ys[isol:isol]
                end
            end
        end
    end
    


    # if kind == (:raoff, :decoff)
    #     xs = raoff.(solns)
    #     ys = decoff.(solns)



    #     xguide --> "Δra - (mas)"
    #     yguide --> "Δdec - (mas)"




    # elseif kind isa Tuple || kind ∈ (:x, :y, :z)
    #     if kind isa Symbol
    #         kind = (kind,)
    #     end
    #     varx = kind[1]
    #     if length(kind) > 1
    #         vary = kind[2]
    #     end
    #     if length(kind) > 2
    #         varz = kind[3]
    #     end
    #     funcs = (;
    #         x=PlanetOrbits.posx,
    #         y=PlanetOrbits.posy,
    #         z=PlanetOrbits.posz,
    #     )
    #     guides = (;
    #         x="x (au)",
    #         y="y (au)",
    #         z="z (au)",
    #     )

    #     xs = funcs[varx].(solns)
    #     x = xs[begin]
    #     if length(kind) > 1
    #         ys = funcs[vary].(solns)
    #         y = ys[begin]
    #     end
    #     if length(kind) > 2
    #         zs = funcs[varz].(solns)
    #         z = zs[begin]
    #     end


    #     # Set the guides if unset
    #     if length(kind) == 1
    #         xguide --> "time (days)"
    #         yguide --> guides[varx]
    #     else
    #         xguide --> guides[varx]
    #         # We almost always want to see spatial coordinates with equal step sizes
    #         aspect_ratio --> 1
    #         xflip --> true
    #     end
    #     if length(kind) > 1
    #         yguide --> guides[vary]
    #     end
    #     if length(kind) > 2
    #         zguide --> guides[varz]
    #     end

    #     if length(kind) == 1
    #         ts = _time_from_soln.(solns)
    #         t = ts[begin]
    #         # Prevent wrapped line
    #         i = findfirst(ts[begin:end-1] .> ts[begin+1:end])
    #         xs[i] = NaN
    #         @series begin
    #             ts, xs
    #         end
    #         if get(plotattributes, :solmarker, true)
    #             @series begin
    #                 seriestype --> :scatter
    #                 label --> ""
    #                 [t], [x]
    #             end
    #         end
    #     elseif length(kind) == 2
    #         @series begin
    #             xs, ys
    #         end
    #         if get(plotattributes, :solmarker, true)
    #             @series begin
    #                 seriestype --> :scatter
    #                 label --> ""
    #                 color --> :gray
    #                 [x], [y]
    #             end
    #         end
    #     elseif length(kind) == 3
    #         @series begin
    #             xs, ys, zs
    #         end
    #         if get(plotattributes, :solmarker, true)
    #             @series begin
    #                 seriestype --> :scatter
    #                 label --> ""
    #                 color --> :gray
    #                 [x],[y],[z]
    #             end
    #         end
    #     else
    #         error("Cannot plot in more than 3 dimensions")
    #     end
    # elseif kind == :radvel
    #     rvs = radvel.(solns)
    #     ts = _time_from_soln.(solns)
    #     # Prevent wrapped line
    #     i = findfirst(ts[begin:end-1] .> ts[begin+1:end])
    #     rvs[i] = NaN

    #     xguide --> "time (days)"
    #     yguide --> "secondary radial velocity (m/s)"

    #     @series begin
    #         ts, rvs
    #     end
    #     if get(plotattributes, :solmarker, true)
    #         @series begin
    #             seriestype --> :scatter
    #             label --> ""
    #             color --> :gray
    #             ts[begin:begin], rvs[begin:begin]
    #         end
    #     end
    # else
    #     error("unrecognized kind=$kind value")
    # end


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


function _time_from_soln(os::AbstractOrbitSolution, mass=nothing; tref=58849)

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

# # https://discourse.julialang.org/t/equivalent-of-matlabs-unwrap/44882/4?
# function unwrap!(x)
# 	v = first(x)
# 	@inbounds for k = eachindex(x)
# 		x[k] = v = v + rem2pi(x[k] - v,  RoundUp)
# 	end
# end
function unwrap!(x, period = 2π)
	y = convert(eltype(x), period)
	v = first(x)
	@inbounds for k = eachindex(x)
		x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
	end
end