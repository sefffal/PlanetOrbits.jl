"""
Module for calculating orbits and converting into cartiesian coordinates.

# Available Functions
- Orbit(...)
- τ2t0(τ,P)
- period(orbit)
- xyz(orbit, t)
- plot(orbit)
"""
module DirectOrbits

using DirectImages


# Required packages
using LinearAlgebra
using CoordinateTransformations
using StaticArrays
using Roots # For solving for eccentric anomaly
import Dates
import Base.inv
using Statistics: mean


# """
# Convert from (fraction of orbit past periastron at MJD=0)
# to the time of periastron passage in MJD.
# P in years.
# τ [0,1]
# """
# function τ2t0(τ,P, τ_ref_epoch=58849) # modern convention as default
#     # t0 = τ_ref_epoch - τ*P*365.25

#     t0 = - τ_ref_epoch - τ
# end



## Function that takes orbital parameters and time, and maps them to an image location
struct Orbit{T<:Real}
    a::T
    i::T
    e::T
    τ::T
    μ::T
    ω::T
    Ω::T
    plx::T
end
Orbit(;a, i, e, τ, μ, ω, Ω, plx) = Orbit(a, i, e, τ, μ, ω, Ω, plx)
export Orbit

import Base: length, iterate
length(::Orbit) = 1
iterate(orb::Orbit) = (orb, nothing)
iterate(::Orbit, ::Nothing) = nothing

"""
Period of orbit in years.
"""
period(orb::Orbit) = sqrt(orb.a^3/orb.μ)
export period

"""
    xyz(orbit, t)
Call an Orbit object with a time (as a modified Julian date) to get
a projected position x,y,and z in milliarcseconds.
"""
function xyz(orb::Orbit, t)
    mas2rad = 4.8481368E-9
    rad2mas = 206264806.2471
    pc2au = 206265e3
    au2m = 1.495978707e11
    
    # Distance to the central object in AU
    dist = 1/(orb.plx*1e-3) * pc2au

    # Compute period (years)
    T = sqrt(orb.a^3/orb.μ)
    n = 2π/T

    # Compute mean anomaly
    # MA = n*(t-orb.t0)/365.25

    # Similar to orbitize!
    # MA = n * (t - 58849)/365.25 - 2π*orb.τ

    # Me messing around
    MA = n * (t - 58849)/365.25 - 2π*orb.τ

    MA = rem2pi(MA, RoundUp)
    
    # Solve for eccentric anomaly
    # The let-block is to capture the bindings of e and M1 directly (performance)
    f, fp = let e = orb.e, MA=MA
        @inline f(EA) = EA - MA - e*sin(EA)
        @inline fp(EA) = 1 - cos(EA)
        f,fp
    end
    # For pathalogical cases, this may not converge.
    # In that case, throw a warning and send the point to the origin (not much else we can do)
    local EA
    try
        # Newton's method
        # EA = Roots.newton(f, fp, π, MA)
        # Bisection method.
        # This seems both slight faster and more robust.
        # EA = find_zero(f, (-π, π))
        EA = find_zero(f, MA)
    catch err
        if typeof(err) <: InterruptException
            rethrow(err)
        end
        @warn "Solving for eccentric anomaly failed" exception=err maxlog=5
        return SVector(0.,0.,0.0)
    end
    # @show EA
    # EA = MA
    # Calculate true anomaly
    ν = 2atan(√((1+orb.e)/(1-orb.e))*tan(EA/2))
    # New orbital radius.
    # This is the semi-major axis, modified by the eccentricity. Units of AO
    r = orb.a*(1-orb.e*cos(EA))
    # h = sqrt(orb.μ*orb.a*(1-orb.e^2))
    # x = r*(cos(orb.Ω)*cos(orb.ω + ν) - sin(orb.Ω)*sin(orb.ω+ν)*cos(orb.i))
    # y = r*(sin(orb.Ω)*cos(orb.ω + ν) - cos(orb.Ω)*sin(orb.ω+ν)*cos(orb.i))
    
    # Project back into Cartesian coordinates (AU).
    x = r*(cos(orb.Ω)*cos(orb.ω+ν) - sin(orb.Ω)*sin(orb.ω+ν)*cos(orb.i))
    y = r*(sin(orb.Ω)*cos(orb.ω+ν) + cos(orb.Ω)*sin(orb.ω+ν)*cos(orb.i))

    # Now convert back to projected separation in pixels
    x = atan(x,dist)*rad2mas # radians -> mas
    y = atan(y,dist)*rad2mas # radians -> mas
    z = atan(r*(sin(orb.i)*sin(orb.ω+ν)),dist)*rad2mas # radians -> mas
    return SVector(y,x,z)#.*1e3
end
export xyz


using RecipesBase
@recipe function f(orb::Orbit)
    ts = range(0, 365.25period(orb), length=200)
    coords = xyz.(orb, ts)
    xs = [c[1] for c in coords].*1e3
    ys = [c[2] for c in coords].*1e3

    # We almost always want to see spatial coordinates with equal step sizes
    aspect_ratio --> 1
    # And we almost always want to reverse the RA coordinate to match how we
    # see it in the sky.
    xflip --> true

    xguide --> "x - mas"
    yguide --> "y - mas"
    framestyle  --> :box
    fontfamily  -->  "Times"
    minorticks  -->  true

    minorgrid --> true
    gridalpha --> 0.4
    return xs, ys
end


end # module
