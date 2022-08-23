#=
This file defines a Coordinate Transform that allows us to project
images forward and backwards in time according to Kepler's laws.
The currently solution is incomplete but workable for inclined
orbital planes.
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
    M::T
    ω::T
    Ω::T
    plx::T

    # Image properties
    platescale::T

    # Time properties
    dt::T


    # Cached constants for these elements.
    dist::T
    ν_fact::T
    cos_Ω::T
    sin_Ω::T
    cos_i::T
    sin_i::T

    # Inner constructor to inforce invariants and pre-calculate a few
    # constants for these elements.
    function OrbitalTransformation(i, e, M, ω, Ω, plx, platescale, dt)


        # Enforce invariants on user parameters
        e = max(zero(e), min(e, one(e)))
        M = max(M, zero(M))
        plx = max(plx, zero(plx))
        # Pre-calculate some factors that will be re-used when calculating kep2cart at any time
        # Distance in AU
        dist = 1/(plx/1000) * pc2au
        # Compute period (days)
        # Factor in calculating the true anomaly
        ν_fact = √((1+e)/(1-e))

        if e != 0
            error("eccentric transformations are not currently working correctly")
        end

        T = promote_type(
            typeof(i), 
            typeof(e),
            typeof(M),
            typeof(ω),
            typeof(Ω),
            typeof(plx),
            typeof(platescale),
            typeof(dt),
        )
        # The user might pass in integers, but it makes no sense to do these
        # calculations on integers. Assume they mean to use floats.
        if T <: Integer
            T = promote_type(T, Float64)
        end

        sin_Ω, cos_Ω = sincos(Ω)
        sin_i, cos_i = sincos(i)
        new{T}(
            # Passed parameters that define the elements
            i,
            e,
            M,
            ω,
            Ω,
            plx,
            platescale,
            dt,
            # Cached calcuations
            dist,            
            ν_fact,
            # Geometric factors
            cos_Ω,
            sin_Ω,
            cos_i,
            sin_i,
        )
    end
end
# Allow arguments to be specified by keyword.
OrbitalTransformation(;i, e, M, ω, Ω, plx, platescale, dt) = OrbitalTransformation(i, e, M, ω, Ω, plx, platescale, dt)
# And by a named tuple without splatting
OrbitalTransformation(nt::NamedTuple) = OrbitalTransformation(nt.i, nt.e, nt.M, nt.ω, nt.Ω, nt.plx, nt.platescale, nt.dt)
export OrbitalTransformation

function (ot::OrbitalTransformation{T})(dist_proj_px) where T
    # Given x and y, solve for a and mean anomaly. Or I guess eccentric anomaly? and then work back?
    # nvm we solve for true anomaly

    # Convert pixel coordinates into AU`
    # dist = 1/(ot.plx/1000) * pc2au

    dist_proj_mas = dist_proj_px*ot.platescale
    dist_proj_as = dist_proj_mas/1e3
    dist_proj_rad = dist_proj_as / rad2as
    dist_proj_au = tan.(dist_proj_rad) .* ot.dist
    (y₀,x₀) = dist_proj_au

    # cos_i = cos(ot.i)
    # sin_Ω, cos_Ω = sincos(ot.Ω)
    # ν_fact = √((1+ot.e)/(1-ot.e))

    r₀′ = √(x₀^2 + y₀^2)

    # Singularity at the origin that has a trivial solution
    if r₀′ ≈ 0
        return SVector{2,typeof(x₀)}(0.0, 0.0)
    end

    # Derive true anomaly from position and known orbital properties
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
    
    
    # ν₀ = atan(
    #     (x₀*sin_Ω - y₀*cos_Ω ),
    #     (y₀*sin_Ω*cos_i - x₀*cos_Ω*cos_i)
    # ) - ω

    # x = r*(sin_Ω*cos_ω_ν + cos_Ω*sin_ω_ν*cos_i)
    # y = r*(cos_Ω*cos_ω_ν - sin_Ω*sin_ω_ν*cos_i)

    # y/x = r*(cos_Ω*cos_ω_ν - sin_Ω*sin_ω_ν*cos_i)/r*(sin_Ω*cos_ω_ν + cos_Ω*sin_ω_ν*cos_i)
    # y/x = (cos_Ω*cos_ω_ν - sin_Ω*sin_ω_ν*cos_i)/
    #       (sin_Ω*cos_ω_ν + cos_Ω*sin_ω_ν*cos_i)
    # y/x = (A*cos_ω_ν - B*sin_ω_ν*C)/ # C = cos_i
    #       (B*cos_ω_ν + A*sin_ω_ν*C)
    # y/x = (A*cos_θ - B*sin_θ*C)/
    #       (B*cos_θ + A*sin_θ*C)
    # y/x = (A*cos_θ - B*sin_θ*C)/
    #       (B*cos_θ + A*sin_θ*C)

    # Substitute:
    # tan θ = sin θ / cos θ
    # sin θ = tan θ * cos θ  

    # y/x = (A*cos_θ - B*tan θ * cos θ *C)/
    #       (B*cos_θ + A*tan θ * cos θ *C)
          
    # y/x = (A - B*tan θ * C)/
    #       (B + A*tan θ * C)

    # y * (B + A*tan θ * C) = x * (A - B*tan θ * C)
    # y * (B + A*Z) = x * (A - B*Z) # Z = tan θ * C
    #     By + A*Zy = Ax - B*Zx
    #    (Ay + Bx)Z = Ax - By
    
    # Expand again
    #           Z = (cos_Ω*x - sin_Ω*y)/(cos_Ω*y + sin_Ω*x)
    #   tan θ * C = (cos_Ω*x - sin_Ω*y)/(cos_Ω*y + sin_Ω*x)
    # tan ω_ν * C = (cos_Ω*x - sin_Ω*y)/(cos_Ω*y + sin_Ω*x)
    #     tan ω_ν = (cos_Ω*x - sin_Ω*y)/(cos_Ω*y + sin_Ω*x)/C
    #     tan ω_ν = (cos_Ω*x - sin_Ω*y)/(cos_Ω*y + sin_Ω*x)  / cos_i
    
    # ω+ν = atan(cos_Ω*x - sin_Ω*y, (cos_Ω*y + sin_Ω*x)/cos_i)
    #   ν = atan(cos_Ω*x - sin_Ω*y, (cos_Ω*y + sin_Ω*x)/cos_i) - ω

    # Code:
    # ν₀ = atan(cos_Ω*x₀ - sin_Ω*y₀, (cos_Ω*y₀ + sin_Ω*x₀)/cos_i) - ot.ω
    
    # De-project
    y₀′ = (ot.cos_Ω*x₀ - ot.sin_Ω*y₀)/ot.cos_i
    x₀′ = (ot.cos_Ω*y₀ + ot.sin_Ω*x₀)

    # Calculate true anomaly of initial position
    ν₀ = atan(y₀′, x₀′) - ot.ω

    # From those, get the initial eccentric anomal
    EA₀ = 2atan(ot.ν_fact*tan(ν₀/2))

    # The inverse of kepler's equation is closed form.
    # This gives us initial mean anomal
    MA₀ = EA₀ - ot.e*sin(EA₀)

    # Orbital radius
    # After we have done the inclination re-projection,
    # x₀′ and y₀′ give the separation from the star in AU
    # and so we can calculate the initial separation as follows
    r₀ = √(x₀′^2 + y₀′^2)

    # Since we now know the eccentric anomaly via the true anomaly,
    # and we have the separation, we can calculate the semi-major axis.
    a = r₀/(1-ot.e*cos(EA₀))

    # Calculate mean motion for this semi-major axis
    m = 2π/√(a^3/ot.M)

    # Advance mean anomaly by dt
    MA = MA₀ + m/convert(T, year2day) * (-1)* ot.dt

    # And finally solve for eccentric anomaly as usual
    EA = kepler_solver(MA, ot.e, Markley())

    ν = convert(T,2)*atan(ot.ν_fact*tan(EA/convert(T,2)))
    
    sin_ω_ν, cos_ω_ν = sincos(ot.ω+ν)
    # sin_Ω, cos_Ω = sincos(ot.Ω)
    r = a*(one(T)-ot.e*cos(EA))
    y = r*(ot.cos_Ω*cos_ω_ν - ot.sin_Ω*sin_ω_ν*ot.cos_i)
    x = r*(ot.sin_Ω*cos_ω_ν + ot.cos_Ω*sin_ω_ν*ot.cos_i)    

    # TODO: move to tests
    if ot.dt == 0
        xtol = isapprox(x,x₀,atol=1e-2)
        ytol = isapprox(y,y₀,atol=1e-2)
        if !xtol || !ytol
            @show x₀ x y₀ y r₀ r a MA MA₀ EA ν₀ ν
            error("Output != input despite dt=0")
        end
    end


    # z = r*(sin(ot.i)*sin(ot.ω+ν))
    coords_AU = SVector(y,x)
    dist_proj_rad = atan.(coords_AU, ot.dist)
    dist_proj_mas = dist_proj_rad .* convert(eltype(dist_proj_rad),rad2as*1e3) # radians -> mas
    dist_proj_px = dist_proj_mas./ot.platescale
    return dist_proj_px

end


# Inverse transform is just time reversal:
Base.inv(orbit::OrbitalTransformation) = OrbitalTransformation(orbit.i,orbit.e,orbit.M,orbit.ω,orbit.Ω,orbit.plx,orbit.platescale,-orbit.dt)