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

@recipe function f(oses::AbstractArray{<:AbstractOrbitSolution})

    label --> ""
    seriesalpha --> 30/length(oses)
    for os in oses
        @series begin
            os
        end
    end
end

default_plotkind(::OrbitSolutionCartesian) = (:x, :y)
default_plotkind(::OrbitSolutionKep) = (:x, :y)
default_plotkind(::OrbitSolutionThieleInnes) = :astrometry
default_plotkind(::OrbitSolutionRadialVelocity) = :radvel
default_plotkind(::OrbitSolutionVisual) = :astrometry


# Plotting recipes for orbital elements
using RecipesBase
@recipe function f(os::AbstractOrbitSolution)

    # Variables to plot
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
    astromaccvars = (:accra, :accdec)
    ravars = (:x, :raoff, :pmra)
    timevars = (:t, :ν, :trueanom, :meananom, :eccanom)
    unwrapvars = (timevars..., :posangle)
    spatial_vars = (posvars..., astromvars..., astromvelvars..., astromaccvars...)
    if all(∈(spatial_vars), kind)
        # Set equal aspect ratio when plotting all spatial coordinates
        aspect_ratio --> 1
    end

    resolver = (;
        t=("t", "mjd", (sol,args...)->_time_from_EA(sol,eccanom(sol))),
        ν=("ν", "rad", trueanom,),
        trueanom=("ν", "rad", trueanom,),
        meananom=("mean.anom.", "rad", meananom,),
        eccanom=("ecc.anom", "rad", eccanom,),
        x=("x", "au", posx,),
        y=("y", "au", posy,), 
        z=("z", "au", posz,),
        velx=("∂x/δt", "au/yr", velx),
        vely=("∂y/δt", "au/yr", vely),
        velz=("∂z/δt", "au/yr", velz),
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
    xguide --> "$xl [$xu]"
    yguide --> "$yl [$yu]"
    if length(kind) >= 3
        zl, zu, zf = resolver[kind[3]]
        zguide --> "$zl [$zu]"
    end

    bodies = get(plotattributes, :body, :secondary)
    if bodies isa Symbol
        bodies = (bodies,)
    end

    L = get(plotattributes, :orbitsteps, 90)
    for body in bodies
        # We trace out in equal steps of eccentric anomaly instead of time for a smooth
        # curve, regardless of eccentricity.
        # When the independent variable is a timevar (angle or time) we want
        # two cycles, otherwise we just do one complete orbit
        if kind[1] == :t || !isfinite(period(os.elem)) 
            if isfinite(period(os.elem))
                # bound orbit case
                default_tspan = (soltime(os)-period(os.elem), soltime(os)+period(os.elem))
                tspan = get(plotattributes, :tspan, default_tspan)
                tstart, tstop = extrema(tspan)
                ea_start = eccanom(orbitsolve(os.elem, tstart))
                ea_stop =  eccanom(orbitsolve(os.elem, tstop))
                while _time_from_EA(os.elem, ea_start) > tstart + 0.000001
                    ea_start -= 2π
                end
                while _time_from_EA(os.elem, ea_stop) < tstop - 0.000001
                    ea_stop += 2π
                end
                # if ea_stop < ea_start
                #     ea_stop += 2π
                # end
                eccanoms = range(
                    ea_start,
                    ea_stop,
                    length=L,
                )
            else
                # non-elliptical case
                default_tspan = (soltime(os)-5*365*meanmotion(os.elem), soltime(os)+5*365*meanmotion(os.elem))
                tspan = get(plotattributes, :tspan, default_tspan)
                tstart, tstop = extrema(tspan)
                eccanoms = range(
                    eccanom(orbitsolve(os.elem, tstart)),
                    eccanom(orbitsolve(os.elem, tstop)),
                    length=L,
                )
            end
        elseif kind[1] ∈ timevars
            eccanoms = range(-2π, 2π, length=L)
            xticks --> (range(-2π, 2π, step=π/2), ["-2π", "-3π/2", "-π", "-π/2", "0", "+π/2", "+π", "+3π/2", "+2π"])
        else
            # Otherwise we are plotting two variables against each other and don't need
            # to consider multiple cycles
            eccanoms = range(eccanom(os), eccanom(os)+2π, length=L)
            line_z --> -eccanoms
            colorbar --> nothing
        end

        if get(plotattributes, :timestep, false)
            if kind[1] == :t
                tspan = get(plotattributes, :tspan, (soltime(os)-period(os.elem), soltime(os)+2period(os.elem)))
                solns = orbitsolve.(os.elem, range(tspan..., length=L))
            else
                solns = orbitsolve.(os.elem, range(0, period(os.elem), length=L))
            end
        else
            solns = orbitsolve_eccanom.(os.elem, eccanoms)
        end

        if body == :secondary
            x = xf(os)
            xs = xf.(solns)
            y = yf(os)
            ys = yf.(solns)
            if length(kind) >= 3
                z = zf(os)
                zs = zf.(solns)
            end
        elseif body == :primary
            x = xf(os, plotattributes[:mass])
            xs = xf.(solns, plotattributes[:mass])
            y = yf(os, plotattributes[:mass])
            ys = yf.(solns, plotattributes[:mass])
            if length(kind) >= 3
                z = zf(os, plotattributes[:mass])
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
        if kind[1] ∈ unwrapvars && isfinite(period(os.elem))
            P = kind[1]==:t ? period(os.elem) : 2π
            unwrap!(xs, P)
            xs .-= P
        end

        @series begin
            label --> string(body)
            if haskey(plotattributes, :seriescolor) 
                line_z := nothing
            end
            if isdefined(Main, :Plots) && isdefined(Main.Plots, :palette) && get(plotattributes, :solmarker, true)
                # We would like to create a nice semi-transparent 
                # gray gradient. But palette isn't in PlotRecipes so
                # we fall back to this hacky way of getting it
                if body == :secondary
                    seriescolor --> Main.Plots.palette(["#444444ff", "#44444433"],10)
                else
                    seriescolor --> Main.Plots.palette(["#BB4444ff", "#BB444433"],10)
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
                    seriescolor --> :gray
                else
                    seriescolor --> "#BB4444"
                end

                if length(kind) >= 3
                    [x], [y], [z]
                else
                    [x], [y]
                end
            end
        end
    end
    

end


# https://discourse.julialang.org/t/equivalent-of-matlabs-unwrap/44882/4?
function unwrap!(x, period = 2π)
	y = convert(eltype(x), period)
	v = first(x)
	@inbounds for k = eachindex(x)
		x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
	end
end
