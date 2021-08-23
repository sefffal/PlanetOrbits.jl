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

    # Convert pixel coordinates into AU`
    dist = 1/(ot.plx/1000) * pc2au
    dist_proj_mas = dist_proj_px*ot.platescale
    dist_proj_as = dist_proj_mas/1e3
    dist_proj_rad = dist_proj_as / rad2as
    dist_proj_au = tan.(dist_proj_rad) .* dist 
    (x₀,y₀) = dist_proj_au

    cos_i = cos(ot.i)
    sin_Ω, cos_Ω = sincos(ot.Ω)
    ν_fact = √((1+ot.e)/(1-ot.e))

    r₀ = √(x₀^2 + y₀^2)

    # Singularity at the origin that has a trivial solution
    if r₀ ≈ 0
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
    # y/x = (A*cos_ω_ν - B*sin_ω_ν*C)/
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
    # y * (B + A*Z) = x * (A - B*Z)
    #     By + A*Zy = Ax - B*Zx
    #    (Ay + Bx)Z = Ax - By
    
    # Expand again
    #           Z = (cos_Ω*x - sin_Ω*y)/(cos_Ω*y + sin_Ω*x)
    #   tan θ * C = (cos_Ω*x - sin_Ω*y)/(cos_Ω*y + sin_Ω*x)
    # tan ω_ν * C = (cos_Ω*x - sin_Ω*y)/(cos_Ω*y + sin_Ω*x)
    #     tan ω_ν = (cos_Ω*x - sin_Ω*y)/(cos_Ω*y + sin_Ω*x)/C
    #     tan ω_ν = (cos_Ω*x - sin_Ω*y)/(cos_Ω*y + sin_Ω*x)/cos_i
    
    # ω+ν = atan(cos_Ω*x - sin_Ω*y, (cos_Ω*y + sin_Ω*x)/cos_i)
    #   ν = atan(cos_Ω*x - sin_Ω*y, (cos_Ω*y + sin_Ω*x)/cos_i) - ω
    ν₀ = atan(cos_Ω*x₀ - sin_Ω*y₀, (cos_Ω*y₀ + sin_Ω*x₀)/cos_i) - ot.ω
          



    EA₀ = 2atan(ν_fact*tan(ν₀/2))

    # We don't use this, but here is the z-coordinate
    # z0 = r0*sin(orb.i)*sin(orb.ω+ν0)

    MA₀ = EA₀ - ot.e*sin(EA₀)
    a = r₀/(1-ot.e*cos(EA₀))

    # Calculate mean motion for this semi-major axis
    m = 2π/√(a^3/ot.μ)

    # Advance mean anomaly by dt
    MA = MA₀ + m/convert(T, year2days) * ot.dt

    # And finally solve as usual
    # MA = rem2pi(MA, RoundDown)
    EA = _kepler_solver_inline(MA, ot.e)

    a = r₀/(1-ot.e*cos(EA))

    ν = convert(T,2)*atan(ν_fact*tan(EA/convert(T,2)))
    

    sin_ω_ν, cos_ω_ν = sincos(ot.ω+ν)
    sin_Ω, cos_Ω = sincos(ot.Ω)
    r = a*(one(T)-ot.e*cos(EA))
    y = r*(cos_Ω*cos_ω_ν - sin_Ω*sin_ω_ν*cos_i)
    x = r*(sin_Ω*cos_ω_ν + cos_Ω*sin_ω_ν*cos_i)

    if ot.dt == 0
        xtol = isapprox(x,x₀,atol=1e-2)
        ytol = isapprox(y,y₀,atol=1e-2)
        if !xtol || !ytol
            @show x₀ x y₀ y r₀ r a MA MA₀ EA ν₀ ν
            error("Output != input despite dt=0")
        end
    end


    # z = r*(sin(ot.i)*sin(ot.ω+ν))
    coords_AU = SVector(x,y)#,z)
    dist_proj_rad = atan.(coords_AU, dist)
    dist_proj_mas = dist_proj_rad .* convert(eltype(dist_proj_rad),rad2as*1e3) # radians -> mas
    dist_proj_px = dist_proj_mas./ot.platescale
    return dist_proj_px

end


# Inverse transform is just time reversal:
Base.inv(orbit::OrbitalTransformation) = OrbitalTransformation(orbit.i,orbit.e,orbit.μ,orbit.ω,orbit.Ω,orbit.plx,orbit.platescale,-orbit.dt)