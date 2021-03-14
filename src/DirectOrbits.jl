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

using LinearAlgebra
# using CoordinateTransformations
using StaticArrays
using Roots # For solving for eccentric anomaly
# import Dates
# import Base.inv
using Statistics: mean


const mas2rad = 4.8481368E-9
const rad2as = 206265
const pc2au = 206265
const au2m = 1.495978707e11
const year2days = 365.2422

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

abstract type AbstractElements end

"""
    Orbit(
        a=1.0, # semi-major axis, AU
        i=π/2, # inclination, radians
        e=0.1, # eccentricity
        τ=π/2, # fraction of orbit past periastron passage at MJD=0,
        μ=1.0, # graviational parameter, solar masses
        ω=π/2, # argument of periapsis
        Ω=π/2, # longitude of the ascending node
        plx=10.1, # paralax in milliarcseconds. Defines the distance to the object
    )

Represents one object's Keplerian orbital elements. Values can be specified
by keyword argument for convinience, or position for speed.
"""
struct KeplarianElements{T<:Number} <: AbstractElements

    # Orbital properties
    a::T
    i::T
    e::T
    τ::T
    μ::T
    ω::T
    Ω::T
    plx::T

    # Cached constants for this orbit
    dist::T
    T::T
    n::T
    ν_fact::T
    cos_Ω::T
    sin_Ω::T
    cos_i::T

    # Inner constructor to inforce invariants, and pre-calculate a few
    # constants for this orbit
    function KeplarianElements(a::T, i::T, e::T, τ::T, μ::T, ω::T, Ω::T, plx::T) where T


        # Enforce invariants on user parameters
        a = max(a, zero(a))
        e = max(zero(e), min(e, one(e)))
        μ = max(μ, zero(μ))
        plx = max(plx, zero(plx))
        # Pre-calculate some factors that will be re-used when calculating position at any time
        # Distance in AU
        dist = 1/(plx*1e-3) * pc2au
        # Compute period (days)
        period = √(a^3/μ)*year2days
        # Mean motion
        n = 2π/√(a^3/μ)
        # Factor in calculating the true anomaly
        ν_fact = √((1+e)/(1-e))
        new{T}(
            # Passed parameters that define the orbit
            a,
            i,
            e,
            τ,
            μ,
            ω,
            Ω,
            plx,
            # Cached calcuations
            dist,            
            period,
            n,
            ν_fact,
            # Geometric factors
            cos(Ω),
            sin(Ω),
            cos(i),
        )
    end
end
# Allow arguments to be specified by keyword.
KeplarianElements(;a, i, e, τ, μ, ω, Ω, plx) = KeplarianElements(a, i, e, τ, μ, ω, Ω, plx)
export KeplarianElements
function Orbit(args...; kwargs...)
    @warn "Orbit is deprecated in favour of KeplerianElements"
    return KeplarianElements(args...; kwrags...)
end
export Orbit

# Better printing
Base.show(io::IO, ::MIME"text/plain", orb::KeplarianElements) = print(
    io, """
        $(typeof(orb))
        ─────────────────────────
        a   [au ] = $(round(orb.a,sigdigits=3)) 
        i   [°  ] = $(round(rad2deg(orb.i),sigdigits=3))
        e         = $(round(orb.e,sigdigits=3))
        τ         = $(round(orb.τ,sigdigits=3))
        μ   [M⊙ ] = $(round(orb.μ,sigdigits=3)) 
        ω   [°  ] = $(round(rad2deg(orb.ω),sigdigits=3))
        Ω   [°  ] = $(round(rad2deg(orb.Ω),sigdigits=3))
        plx [mas] = $(round(orb.plx,sigdigits=3)) 
        ──────────────────────────
        period      [yrs ] : $(round(period(orb)/year2days,digits=1)) 
        distance    [pc  ] : $(round(distance(orb),digits=1)) 
        mean motion [°/yr] : $(round(rad2deg(meanmotion(orb)),sigdigits=3)) 
        ──────────────────────────
        """)
Base.show(io::IO, orb::KeplarianElements) = print(io,
    "KeplarianElements($(round(orb.a,sigdigits=3)), $(round(orb.i,sigdigits=3)), $(round(orb.e,sigdigits=3)), "*
    "$(round(orb.τ,sigdigits=3)), $(round(orb.μ,sigdigits=3)), $(round(orb.ω,sigdigits=3)), "*
    "$(round(orb.Ω,sigdigits=3)), $(round(orb.plx,sigdigits=3)))"
)

import Base: length, iterate
length(::AbstractElements) = 1
iterate(orb::AbstractElements) = (orb, nothing)
iterate(::AbstractElements, ::Nothing) = nothing

"""
Period of orbit in days.
"""
period(orb::KeplarianElements) = orb.T
export period

"""
Distance to the system in parsecs.
"""
distance(orb::KeplarianElements) = orb.dist/pc2au
export distance

"""
Mean motion, radians per year.
"""
meanmotion(orb::KeplarianElements) = orb.n

"""
    xyz(orbit, t, [throw_ea=false])

Given an Orbit value with a time (as a modified Julian day) to get
a projected position x, y, and z in milliarcseconds.

In pathalogical cases solving for eccentric anomaly might fail.
This is unlikely for any reasonable orbit, but if using this routine
as part of an image distortion step (via e.g. CoordinateTransformations)
than this can occur near the origin. A warning will be generated
and the function will return (0,0,0). Specifying `throw_ea=true`
turns that warning into an error.
"""
function xyz(orb::KeplarianElements{T}, t, throw_ea=false) where T

    


    # Compute mean anomaly
    # MA = orb.n * (t - 58849.0)/year2days - 2π*orb.τ
    # MA = orb.n * (t/365.25) - 2π*orb.τ

    # MA = orb.τ - orb.n * (t - 58849.0)/year2days 
    # MA = orb.τ - orb.n * t /year2days 
    # MA = orb.n * (t - 58849.0)/365.25 - 2π*orb.τ

    # MA = 2π*orb.τ - orb.n * (t - 58849.0)/year2days 
    # MA = orb.n * (58849.0 - t)/year2days - 2π*orb.τ

    MA = meanmotion(orb)/year2days * (t - orb.τ)

    MA = rem2pi(MA, RoundDown)

    EA = eccentric_anomaly(orb.e, MA; throw_ea)
    
    
    # @show EA
    # EA = MA
    # Calculate true anomaly
    ν = 2atan(orb.ν_fact*tan(EA/2))

    # New orbital radius.
    # This is the semi-major axis, modified by the eccentricity. Units of AO
    r = orb.a*(1-orb.e*cos(EA))


    ## From oribitize for testing


    # # compute ra/dec offsets (size: n_orbs x n_dates)
    # # math from James Graham. Lots of trig
    # c2i2 = cos(0.5*orb.i)^2
    # s2i2 = sin(0.5*orb.i)^2
    # arg1 = ν + orb.ω + orb.Ω
    # arg2 = ν + orb.ω - orb.Ω
    # c1 = cos(arg1)
    # c2 = cos(arg2)
    # s1 = sin(arg1)
    # s2 = sin(arg2)

    # # updated sign convention for Green Eq. 19.4-19.7
    # raoff = r * (c2i2*s1 - s2i2*s2) * orb.plx
    # deoff = r * (c2i2*c1 + s2i2*c2) * orb.plx
    # return SVector(raoff,deoff,0) # as -> mas

    # ##



    # h = sqrt(orb.μ*orb.a*(1-orb.e^2))
    
    # Project back into Cartesian coordinates (AU).
    x = r*(orb.cos_Ω*cos(orb.ω+ν) - orb.sin_Ω*sin(orb.ω+ν)*orb.cos_i)
    y = r*(orb.sin_Ω*cos(orb.ω+ν) + orb.cos_Ω*sin(orb.ω+ν)*orb.cos_i)

    # Now convert back to projected separation in pixels
    x = atan(x,orb.dist)*rad2as # radians -> as
    y = atan(y,orb.dist)*rad2as # radians -> as
    z = atan(r*(sin(orb.i)*sin(orb.ω+ν)),orb.dist)*rad2as # radians -> ma
    
    
    return SVector(y*1e3,x*1e3,z*1e3) # as -> mas
end
export xyz


"""
    eccentric_anomaly(orb, MA; throw_ea)

From an orbit and mean anomaly, calculate the eccentric anomaly
numerically.

In pathalogical cases solving for eccentric anomaly might fail.
This is unlikely for any reasonable orbit, but if using this routine
as part of an image distortion step (via e.g. CoordinateTransformations)
than this can occur near the origin. A warning will be generated
and the function will return (0,0,0). Specifying `throw_ea=true`
turns that warning into an error.
"""
function eccentric_anomaly(e, MA; throw_ea)

    # Numerically solve EA = MA + e * sin(EA) for EA, given MA and e.



    # Fast path for perfectly circular orbits
    if e == 0.0
        return MA
    end

    # Solve for eccentric anomaly
    # The let-block is to capture the bindings of e and M1 directly (performance)
    f= let e = e, MA=MA
        @inline f(EA) = EA - MA - e*sin(EA)
    end

    # After experimentation, Roots finds the root the 
    # fastest / with least allocations using zeroth-order 
    # methods without derivatives. This was surprising.
    # Therefore, we only pass the ojective and no
    # derivatives even though they are trivial to provide.

    # For pathalogical cases, this may not converge.
    # In that case, throw a warning and send the point to the origin


    # For cases very close to one, use a method based on the bisection
    # method
    if isapprox(e, 1, rtol=1e-4)
        try
            # This is a modification of the bisection method. It should be 
            # very very robust.
            EA = find_zero(f, (-2π, 0), FalsePosition(), maxevals=100)
        catch err
            if typeof(err) <: InterruptException
                rethrow(err)
            end
            if throw_ea
                error("Solving for eccentric anomaly near 1 failed")
            else
                @warn "Solving for eccentric anomaly near 1 failed. Pass `throw_ea=true` to turn this into an error." exception=err maxlog=5
                return MA
            end
        end
    end

    # In general, on the other hand:
    local EA
    try
        EA₀ = MA
        # Begin the initial conditions differntly for highly eccentric orbits
        if e > 0.8
            EA₀ = π
        end
        EA = find_zero(f, EA₀, maxevals=100)
    catch err
        if typeof(err) <: InterruptException
            rethrow(err)
        end
        # If it fails to converge in some pathalogical case,
        # try a different root finding algorithm that may work
        try
            # This is a modification of the bisection method. It should be 
            # very very robust.
            EA = find_zero(f, (-2π, 0), FalsePosition(), maxevals=100)
        catch err
            if typeof(err) <: InterruptException
                rethrow(err)
            end
            if throw_ea
                error("Solving for eccentric anomaly failed twice")
            else
                @warn "Solving for eccentric anomaly failed twice. Pass `throw_ea=true` to turn this into an error." exception=err maxlog=5
                return MA
            end
        end
    end
end


# using ChainRulesCore
# # @scalar_rule eccentric_anomaly(y, x) @setup(u = x ^ 2 + y ^ 2) (x / u, -y / u)
# @scalar_rule eccentric_anomaly(e, MA) @setup(u = x ^ 2 + y ^ 2) (x / u, -y / u)

# # function frule((Δself, Δargs...), ::typeof(foo), args...; kwargs...)
# #     ...
# #     return y, ∂Y
# # end


"""
    x(orbit, t)
Call an AbstractElements object with a time (as a modified Julian date) to get
a projected position x.
"""
function x(orb::AbstractElements, t)
    return xyz(orb::AbstractElements, t)[1]
end


"""
    y(orbit, t)
Call an AbstractElements object with a time (as a modified Julian date) to get
a projected position y.
"""
function y(orb::AbstractElements, t)
    return xyz(orb::AbstractElements, t)[2]
end


"""
    z(orbit, t)
Call an AbstractElements object with a time (as a modified Julian date) to get
a projected position z.
"""
function z(orb::AbstractElements, t)
    return xyz(orb::AbstractElements, t)[3]
end

"""
    projectedseparation(orbit, t)

Projected separation in mas from the star at time t (mjd).
"""
function projectedseparation(orb::AbstractElements, t)
    x,y,z = xyz(orb,t)
    return sqrt(x^2 + y^2 + z^2)
end
export projectedseparation


using RecipesBase
@recipe function f(orb::AbstractElements)
    ts = range(0, period(orb), step=year2days/12)
    if length(ts) < 60
        ts = range(0, period(orb), length=60)
    end
    coords = xyz.(orb, ts)
    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]

    # We almost always want to see spatial coordinates with equal step sizes
    aspect_ratio --> 1
    # And we almost always want to reverse the RA coordinate to match how we
    # see it in the sky.
    xflip --> true


    return xs, ys
end

@recipe function f(orbs::AbstractArray{<:AbstractElements})
    ts = range(0, maximum(period.(orbs)), step=year2days/12)
    if length(ts) < 60
        ts = range(0, maximum(period.(orbs)), length=60)
    end
    coords = xyz.(orbs, ts')
    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]

    # We almost always want to see spatial coordinates with equal step sizes
    aspect_ratio --> 1
    # And we almost always want to reverse the RA coordinate to match how we
    # see it in the sky.
    xflip --> true
    xguide --> "ΔRA - mas"
    yguide --> "ΔDEC - mas"

    seriesalpha --> 10/length(orbs)

    return xs, ys
end

# include("Fitting.jl")

end # module
