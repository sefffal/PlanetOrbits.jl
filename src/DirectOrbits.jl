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


const mas2rad = 4.8481368E-9
const rad2as = 206265
const pc2au = 206265
const au2m = 1.495978707e11


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

abstract type AbstractOrbit end

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
struct Orbit{T<:AbstractFloat} <: AbstractOrbit

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
    function Orbit(a::T, i::T, e::T, τ::T, μ::T, ω::T, Ω::T, plx::T) where T <: AbstractFloat
        # Enforce invariants on user parameters
        a = max(a, zero(a))
        e = max(zero(e), min(e, one(e)))
        μ = max(μ, zero(μ))
        plx = max(plx, zero(plx))
        # T = typeof(a*i*e*τ*μ*ω*Ω*plx)
        factor1 = 1/(plx*1e-3) * pc2au
        factor2 = sqrt(a^3/μ)
        factor3 = 2π/sqrt(a^3/μ)
        factor4 = √((1+e)/(1-e))
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
            # Cache constants for this orbit to save time in calculations
            # Distance in AU
            factor1,
            # Compute period (years)
            factor2,
            # Mean motion
            factor3,
            # Factor in calculating the true anomaly
            factor4,
            # Geometric factors
            cos(Ω),
            sin(Ω),
            cos(i),
        )
    end
end
# Allow arguments to be specified by keyword.
Orbit(;a, i, e, τ, μ, ω, Ω, plx) = Orbit(a, i, e, τ, μ, ω, Ω, plx)
export Orbit

# Better printing
Base.show(io::IO, ::MIME"text/plain", orb::Orbit) = print(
    io, """
        Orbit{$(eltype(orb.a))}
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
        period      [yrs ] : $(round(orb.T,digits=1)) 
        distance    [pc  ] : $(round(orb.dist/pc2au,digits=1)) 
        mean motion [°/yr] : $(round(rad2deg(orb.n),sigdigits=3)) 
        ──────────────────────────
        """)
Base.show(io::IO, orb::Orbit) = print(io,
    "Orbit($(round(orb.a,sigdigits=3)), $(round(orb.i,sigdigits=3)), $(round(orb.e,sigdigits=3)), "*
    "$(round(orb.τ,sigdigits=3)), $(round(orb.μ,sigdigits=3)), $(round(orb.ω,sigdigits=3)), "*
    "$(round(orb.Ω,sigdigits=3)), $(round(orb.plx,sigdigits=3)))"
)

import Base: length, iterate
length(::Orbit) = 1
iterate(orb::Orbit) = (orb, nothing)
iterate(::Orbit, ::Nothing) = nothing

"""
Period of orbit in years.
"""
# period(orb::Orbit) = sqrt(orb.a^3/orb.μ)
period(orb::Orbit) = orb.T
export period

"""
    xyz(orbit, t)

Given an Orbit value with and a time (as a modified Julian day) to get
a projected position x, y, and z in milliarcseconds.
"""
function xyz(orb::Orbit, t)

    


    # Compute mean anomaly
    # MA = orb.n * (t - 58849.0)/365.2422 - 2π*orb.τ
    # MA = orb.n * (t/365.25) - 2π*orb.τ

    # MA = orb.τ - orb.n * (t - 58849.0)/365.2422 
    # MA = orb.τ - orb.n * t /365.2422 
    # MA = orb.n * (t - 58849.0)/365.25 - 2π*orb.τ

    # MA = 2π*orb.τ - orb.n * (t - 58849.0)/365.2422 
    MA = orb.n * (58849.0 - t)/365.2422 - 2π*orb.τ

    MA = rem2pi(MA, RoundDown)
    
    # Solve for eccentric anomaly
    # The let-block is to capture the bindings of e and M1 directly (performance)
    f= let e = orb.e, MA=MA
        @inline f(EA) = EA - MA - e*sin(EA)
    end
    fp(EA) = 1 - cos(EA)
    fpp = sin
    # f = let e = orb.e,
    #         MA=MA
    #     @inline f(EA) = EA - MA - e*sin(EA)
    # end
    
    # For pathalogical cases, this may not converge.
    # In that case, throw a warning and send the point to the origin (not much else we can do)
    local EA
    try
        EA = find_zero(f, MA, maxevals=100)
    catch err
        if typeof(err) <: InterruptException
            rethrow(err)
        end
        try
            EA = find_zero(f, (-2π, 0), FalsePosition(), maxevals=100)
        catch err
            if typeof(err) <: InterruptException
                rethrow(err)
            end
            @warn "Solving for eccentric anomaly failed twice" exception=err maxlog=5
            return SVector(0.,0.,0.0)
        end
    end
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
    
    
    return SVector(y,x,z).*1e3 # as -> mas


end
export xyz


"""
    x(orbit, t)
Call an Orbit object with a time (as a modified Julian date) to get
a projected position x.
"""
function x(orb::Orbit, t)
    return xyz(orb::Orbit, t)[1]
end


"""
    y(orbit, t)
Call an Orbit object with a time (as a modified Julian date) to get
a projected position x.
"""
function y(orb::Orbit, t)
    return xyz(orb::Orbit, t)[2]
end


using RecipesBase
@recipe function f(orb::Orbit)
    ts = range(0, 365.25period(orb), step=365.25/4)
    coords = xyz.(orb, ts)
    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]

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

@recipe function f(orbs::AbstractArray{<:Orbit})
    ts = range(0, 365.25maximum(period.(orbs)), step=365.25/4)
    coords = xyz.(orbs, ts')
    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]

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

    seriesalpha --> 10/length(orbs)

    return xs, ys
end

include("Fitting.jl")

end # module
