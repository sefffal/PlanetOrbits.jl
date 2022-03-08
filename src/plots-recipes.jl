#=
This file contains plot recipes for the Plots.jl
ecosystem. This way you can do e.g.:
plot(elems)
=#

# Plotting recipes for orbital elementst
using RecipesBase
@recipe function f(elem::AbstractOrbit)
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
@recipe function f(elems::AbstractArray{<:AbstractOrbit})

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
@recipe function f(os::OrbitSolution)

    color --> :gray
    c = plotattributes[:color]

    # Find true anomaly of orbit solution. There's probably an easier way by 
    # finding the angles of the os coordinates.
        # # Epoch of periastron passage
        # tₚ = periastron(elem, tref)

        # # Mean anomaly    
        # MA = meanmotion(elem)/convert(T2, year2day) * (t - tₚ)
    
        # # Compute eccentric anomaly
        # EA = kepler_solver(MA, elem.e)
        
        # # Calculate true anomaly
        # ν = convert(T2,2)*atan(elem.ν_fact*tan(EA/convert(T2,2)))
    ν = atan(os.x, os.y)+π

    @series begin
        # Hacky
        if isdefined(Main, :Plots)
            color := Main.Plots.palette(["#444444ff", "#44444433"],10)
        end
        
        # We trace out in equal steps of true anomaly instead of time for a smooth
        # curve, regardless of eccentricity.
        νs = range(-π+ν, π+ν, length=90)
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
        color := c
        seriestype --> :scatter
        label --> ""

        [raoff(os)], [decoff(os)]
    end
end