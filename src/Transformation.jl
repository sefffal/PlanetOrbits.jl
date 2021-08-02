#=
Simpler method: just rotate images and stack according
to a given orbit.

INPUT: method, elements, images, epochs
OUTPUT: stacked image
=#


#=
Goals:

I want to be able to enter orbital parameters and get a coordinate transformations
object. 
I guess this object must encode the delta T and the image platescale.
Main use case is specifying everything except semi-major axis.

Then, I want a convenience function for taking an image, orbital parameters, and dt 
and returning a new image.

Finally, I want convenience functions for doing this to a sequence of images and 
returning a stack (custom stacking function of course). Stacking function must
look back into contrast map from original image. Reverse the dt at current location
and lookup via projected separation.

This is *sort* of inverse to what we have now. We need to lookup pixels, and see 
where those pixels were before.

Finally, I want to wrap this up with some animation functions. We should be
able to input an image and step through dt.
Or input multiple images and fade between them.

=#


struct OrbitalTransformation{T<:Number}
    # Orbital properties
    i::T
    e::T
    μ::T
    ω::T
    Ω::T
    plx::T

    # Image properties
    platescale::T

    # Time properties
    dt::T
end
# Allow arguments to be specified by keyword.
OrbitalTransformation(;i, e, μ, ω, Ω, plx, platescale, dt) = OrbitalTransformation(i, e, μ, ω, Ω, plx, platescale, dt)
# And by a named tuple without splatting
OrbitalTransformation(nt::NamedTuple) = OrbitalTransformation(nt.i, nt.e, nt.μ, nt.ω, nt.Ω, nt.plx, nt.platescale, nt.dt)
export OrbitalTransformation

function (ot::OrbitalTransformation{T})(dist_proj_px) where T
    # Given x and y, solve for a and mean anomaly. Or I guess eccentric anomaly? and then work back?
    # nvm we solve for true anomaly

    # Algebra to figure out the right equation for this transformation
    # x = r*(elem.sin_Ω*cos(elem.ω+ν) + elem.cos_Ω*sin(elem.ω+ν)*elem.cos_i)
    # y = r*(elem.cos_Ω*cos(elem.ω+ν) - elem.sin_Ω*sin(elem.ω+ν)*elem.cos_i)
    # x/r = elem.sin_Ω*cos(elem.ω+ν) + elem.cos_Ω*sin(elem.ω+ν)*elem.cos_i
    # y/r = elem.cos_Ω*cos(elem.ω+ν) - elem.sin_Ω*sin(elem.ω+ν)*elem.cos_i
    # y/r + elem.sin_Ω*sin(elem.ω+ν)*elem.cos_i = elem.cos_Ω*cos(elem.ω+ν)
    # (y/r + elem.sin_Ω*sin(elem.ω+ν)*elem.cos_i)/elem.cos_Ω = cos(elem.ω+ν)
    # x/r = elem.sin_Ω*(y/r + elem.sin_Ω*sin(elem.ω+ν)*elem.cos_i)/elem.cos_Ω + elem.cos_Ω*sin(elem.ω+ν)*elem.cos_i
    # x/r*elem.cos_Ω = elem.sin_Ω*(y/r + elem.sin_Ω*sin(elem.ω+ν)*elem.cos_i) + elem.cos_Ω*sin(elem.ω+ν)*elem.cos_i*elem.cos_Ω
    # x/r*elem.cos_Ω/elem.sin_Ω = (y/r + elem.sin_Ω*sin(elem.ω+ν)*elem.cos_i) + elem.cos_Ω*sin(elem.ω+ν)*elem.cos_i*elem.cos_Ω/elem.sin_Ω
    # x/r*elem.cos_Ω/elem.sin_Ω - y/r = elem.sin_Ω*sin(elem.ω+ν)*elem.cos_i + elem.cos_Ω*sin(elem.ω+ν)*elem.cos_i*elem.cos_Ω/elem.sin_Ω
    # x/r*elem.cos_Ω/elem.sin_Ω - y/r = sin(elem.ω+ν)*(elem.sin_Ω*elem.cos_i + elem.cos_Ω*elem.cos_i*elem.cos_Ω/elem.sin_Ω)
    # (x/r*elem.cos_Ω/elem.sin_Ω - y/r)/(elem.sin_Ω*elem.cos_i + elem.cos_Ω*elem.cos_i*elem.cos_Ω/elem.sin_Ω) = sin(elem.ω+ν)
    # (x/r*elem.cos_Ω/elem.sin_Ω - y/r)/(elem.sin_Ω*elem.cos_i + elem.cos_Ω*elem.cos_i*elem.cos_Ω/elem.sin_Ω) = 
    # elem.ω+ν = asin(x/r*elem.cos_Ω/elem.sin_Ω - y/r)/(elem.sin_Ω*elem.cos_i + elem.cos_Ω*elem.cos_i*elem.cos_Ω/elem.sin_Ω)
    # ν = asin(x/r*elem.cos_Ω/elem.sin_Ω - y/r)/(elem.sin_Ω*elem.cos_i + elem.cos_Ω*elem.cos_i*elem.cos_Ω/elem.sin_Ω) - elem.ω

    dist = 1/(ot.plx/1000) * pc2au

    dist_proj_mas = dist_proj_px*ot.platescale
    dist_proj_as = dist_proj_mas/1e3

    # coords_AU = SVector(x,y,z)
    # # coords_AU = MVector(x,y,z)
    # # coords_AU = [x,y,z]
    # dist_proj_rad = atan.(coords_AU, elem.dist)
    # dist_proj_mas = dist_proj_rad .* convert(eltype(dist_proj_rad),rad2as*1e3) # radians -> mas

    # (x₀,y₀) = dist_proj_mas ./ (rad2as*1e3)
    dist_proj_rad = dist_proj_as / rad2as
    dist_proj_au = tan.(dist_proj_rad) .* dist 
    (x₀,y₀) = dist_proj_au

    # @show dist_proj_px dist_proj_as dist_proj_rad dist_proj_au
    
    cos_i = cos(ot.i)
    cos_Ω = cos(ot.Ω)
    sin_Ω = sin(ot.Ω)
    ν_fact = √((1+ot.e)/(1-ot.e))
    ν_fact⁻¹ = √((1-ot.e)/(1+ot.e))

    r₀ = √(x₀^2 + y₀^2)

    @show r₀
    # Singularity at the origin that has a trivial solution
    if r₀ ≈ 0
        return SVector{typeof(x₀)}(0.0, 0.0)
    end

    # # Avoid singularity at cos_Ω==0
    # sin_Ω_safe = cos_Ω ≈ 0 ? one(sin_Ω) : sin_Ω
    # @show x₀ y₀ r₀ cos_i cos_Ω sin_Ω #sin_Ω_safe

    # ν₀ = asin(
    #     (x₀/r₀*cos_Ω/sin_Ω - y₀/r₀)/
    #         (sin_Ω*cos_i + cos_Ω*cos_i*cos_Ω/sin_Ω)
    # ) - ot.ω

    # c1 = -(x₀*sin(ot.Ω)-y₀*cos(ot.Ω))
    # c2 = x₀*cos(ot.Ω)*cos(ot.i) + y₀*sin(ot.Ω)*cos(ot.i)
    # ν₀ = atan(c2,c1) - ot.ω

    
    # c1 = -(y₀*sin(ot.Ω)-x₀*cos(ot.Ω))
    # c2 = y₀*cos(ot.Ω)*cos(ot.i) + x₀*sin(ot.Ω)*cos(ot.i)
    # ν₀ = atan(c1,c2) - ot.ω

    Ω = ot.Ω
    ω = ot.ω
    i = ot.i



    
    ## New approach
    # y/x =  (sin_Ω*cos(θ) + cos_Ω*cos_i*sin(θ)) / 
    #        (cos_Ω*cos(θ) + sin_Ω*cos_i*sin(θ))
    # y/x = (sin_Ω + cos_Ω*cos_i*tan(θ)) / 
    #       (cos_Ω + sin_Ω*cos_i*tan(θ))
    # y*(cos_Ω + sin_Ω*cos_i*tan(θ)) = x*(sin_Ω + cos_Ω*cos_i*tan(θ))    
    # y*cos_Ω + y*sin_Ω*cos_i*tan(θ) = x*sin_Ω + x*cos_Ω*cos_i*tan(θ)
    # y*sin_Ω*cos_i*tan(θ) - x*cos_Ω*cos_i*tan(θ) = x*sin_Ω - y*cos_Ω 
    # tan(θ)*(y*sin_Ω*cos_i - x*cos_Ω*cos_i) = x*sin_Ω - y*cos_Ω 
    # tan(θ) = (x*sin_Ω - y*cos_Ω)/(y*sin_Ω*cos_i - x*cos_Ω*cos_i)
    # ω + ν = atan((x*sin_Ω - y*cos_Ω)/(y*sin_Ω*cos_i - x*cos_Ω*cos_i))
    # ν₀ = atan(x₀*sin_Ω - y₀*cos_Ω, y₀*sin_Ω*cos_i - x₀*cos_Ω*cos_i) - ω
    ν₀ = atan((x₀*sin_Ω - y₀*cos_Ω )/( y₀*sin_Ω*cos_i - x₀*cos_Ω*cos_i)) - ω




    # Semi-working approach.
    # But r is clearly not right...
    # ν₀ = acos((x₀/r₀ - y₀/r₀)/(sin(Ω)*cos(i) + cos(Ω)*cos(i) + sin(Ω) + cos(Ω))) - ω
    # ν₀ = acos(2(x₀/r₀ - y₀/r₀)/(sin(Ω)*cos(i) + cos(Ω)*cos(i) + sin(Ω) + cos(Ω)))/2 - ω




    # # Handle mapping in different quadrants
    # if sign(x₀) == -1 && sign(y₀) == -1
    # if sign(x₀) == -1 && (sign(y₀) == -1 || y₀ ≈ 0)
    if (sign(x₀) == -1 || x₀ ≈ 0) && (sign(y₀) == -1 || y₀ ≈ 0)
        # ν₀ = -π + ν₀ 
    elseif (sign(y₀) == -1 || y₀ ≈ 0)
        # ν₀ = π - ν₀ 
        # ν₀ = 0.
    end

    # # Now we have true anomaly, convert back to mean anomaly.
    # # Then we can go about solving as usual
    # # ν = 2*atan(elem.ν_fact*tan(EA/2))
    # # ν/2 = atan(elem.ν_fact*tan(EA/2))
    # # tan(ν/2) = elem.ν_fact*tan(EA/2)
    # # tan(ν/2)/elem.ν_fact = tan(EA/2)
    # # atan(tan(ν/2)/elem.ν_fact) = EA/2
    # # 2atan(tan(ν/2)/elem.ν_fact) = EA
    # # 2atan(tan(ν/2)/elem.ν_fact⁻¹) = EA
    EA₀ = 2atan(ν_fact*tan(ν₀/2))


    radius = y₀/(sin_Ω*cos(ot.ω+ν₀) + cos_Ω*sin(ot.ω+ν₀)*cos_i)
    # This blows up at y=0, so switch to the alternative formula using x:
    # if y₀ ≈ 0.
    if y₀ < x₀
        radius  = x₀/(cos_Ω*cos(ot.ω+ν₀) - sin_Ω*sin(ot.ω+ν₀)*cos_i)
    end

    # We don't use this, but here is the z-coordinate
    # z0 = r0*sin(orb.i)*sin(orb.ω+ν0)

    # I think the problem might be here.
    # Finally, use these to get the semi-major axis.
    # a = r0*(1+orb.e*cos(ν0))/(1+orb.e^2)
    # a = radius/(1-orb.e*cos(E0))
    a = radius*(1+ot.e*cos(ν₀))/(1+ot.e^2)

    @show a

    # @show EA₀ 

    # Now calculate mean anomaly from eccentric anomaly
    # 0 = EA - MA - e*sin(EA)
    MA₀ = EA₀ - ot.e*sin(EA₀)

    # @show MA₀ 

    # We need to calculate semi-major axis
    # r = elem.a*(1-e*cos(EA))

    # a = r/(1-ecos(EA))

    # r/(1-e*cos(EA)) = a
    # a = r/(1-e*cos(EA)) 
    a = r₀/(1-ot.e*cos(EA₀))

    # @show a

    # Calculate mean motion for this semi-major axis
    m = 2π/√(a^3/ot.μ)

    # @show m

    # Advance mean anomaly by dt
    MA = MA₀ + m/convert(T, year2days) * ot.dt

    # @show m

    # And finally solve as usual
    MA = rem2pi(MA, RoundDown)
    EA = kepler_solver(MA, ot.e)
    ν = convert(T,2)*atan(ν_fact*tan(EA/convert(T,2)))
    r = a*(one(T)-ot.e*cos(EA))
    x = r*(sin_Ω*cos(ot.ω+ν) + cos_Ω*sin(ot.ω+ν)*cos_i)
    y = r*(cos_Ω*cos(ot.ω+ν) - sin_Ω*sin(ot.ω+ν)*cos_i)

    # z = r*(sin(ot.i)*sin(ot.ω+ν))
    coords_AU = SVector(x,y)#,z)
    dist_proj_rad = atan.(coords_AU, dist)
    dist_proj_mas = dist_proj_rad .* convert(eltype(dist_proj_rad),rad2as*1e3) # radians -> mas
    dist_proj_px = dist_proj_mas./ot.platescale
    return dist_proj_px

end


# Inverse transform is just time reversal:
Base.inv(orbit::OrbitalTransformation) = OrbitalTransformation(i,e,μ,ω,Ω,plx,platescale,-dt)