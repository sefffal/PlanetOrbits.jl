"""
# DirectOrbits
A package for calculating orbits in the context of direct imaging.

"""
module DirectOrbits

using LinearAlgebra
# using CoordinateTransformations
using StaticArrays
# using Roots # For solving for eccentric anomaly
# import Dates
# import Base.inv

using AstroLib: kepler_solver


const mas2rad = 4.8481368E-9
const rad2as = 206265
const pc2au = 206265
const au2m = 1.495978707e11
const year2days = 365.2422


using ComponentArrays

abstract type AbstractElements end



"""
    Orbit(
        a=1.0, # semi-major axis, AU
        i=π/2, # inclination, radians
        e=0.1, # eccentricity
        τ=π/2, # fraction of elements past periastron passage at MJD=0,
        μ=1.0, # graviational parameter, solar masses
        ω=π/2, # argument of periapsis
        Ω=π/2, # longitude of the ascending node
        plx=10.1, # paralax in milliarcseconds. Defines the distance to the object
    )

Represents one object's Keplerian elementsal elements. Values can be specified
by keyword argument for convinience, or kep2cart for efficiency.

See also `KeplerianElementsDeg` for a convinience constructor accepting
units of degrees instead of radians.
"""
struct KeplerianElements{T<:Number} <: AbstractElements

    # Orbital properties
    a::T
    i::T
    e::T
    τ::T
    μ::T
    ω::T
    Ω::T
    plx::T

    # Cached constants for these elements.
    dist::T
    T::T
    n::T
    ν_fact::T
    cos_Ω::T
    sin_Ω::T
    cos_i::T
    sin_i::T

    # Inner constructor to inforce invariants and pre-calculate a few
    # constants for these elements.
    function KeplerianElements(a, i, e, τ, μ, ω, Ω, plx)


        if a < 0.0
            @warn "Invalid semi-major axis" a maxlog=50
        end

        if !(0 <= e < 1)
            @warn "Eccentricity out of range" e maxlog=50
        end
        if μ < 0.0
            @warn "Invalid primary mass (<0.001 Msun)" μ maxlog=50
        end

        # Enforce invariants on user parameters
        a = max(a, zero(a))
        e = max(zero(e), min(e, one(e)))
        μ = max(μ, zero(μ))
        plx = max(plx, zero(plx))
        # Pre-calculate some factors that will be re-used when calculating kep2cart at any time
        # Distance in AU
        dist = 1/(plx/1000) * pc2au
        # Compute period (days)
        period = √(a^3/μ) * year2days
        # Mean motion
        n = 2π/√(a^3/μ)
        # Factor in calculating the true anomaly
        ν_fact = √((1+e)/(1-e))

        τ = mod(τ, one(τ))

        T = promote_type(
            typeof(a),
            typeof(i), 
            typeof(e),
            typeof(τ),
            typeof(μ),
            typeof(ω),
            typeof(Ω),
            typeof(plx),
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
            cos_Ω,
            sin_Ω,
            cos_i,
            sin_i,
        )
    end
end
# Allow arguments to be specified by keyword.
KeplerianElements(;a, i, e, τ, μ, ω, Ω, plx) = KeplerianElements(a, i, e, τ, μ, ω, Ω, plx)
# And by a named tuple without splatting
KeplerianElements(nt) = KeplerianElements(nt.a, nt.i, nt.e, nt.τ, nt.μ, nt.ω, nt.Ω, nt.plx)

export KeplerianElements

"""
    astuple(elements)

Return the parameters of a KeplerianElements value as a tuple.
"""
function astuple(elem::KeplerianElements)
    return (;elem.a,elem.i,elem.e,elem.τ,elem.μ,elem.ω,elem.Ω,elem.plx)
end

"""
    KeplerianElementsDeg(a, i, e, τ, μ, ω, Ω, plx)

A convinience function for constructing KeplerianElements where
`i`, `ω`, and `Ω` are provided in units of degrees instead of radians.
"""
KeplerianElementsDeg(a, i, e, τ, μ, ω, Ω, plx) = KeplerianElements(a, deg2rad(i), e, τ, μ, deg2rad(ω), deg2rad(Ω), plx)
KeplerianElementsDeg(;a, i, e, τ, μ, ω, Ω, plx) = KeplerianElementsDeg(a, i, e, τ, μ, ω, Ω, plx)
export KeplerianElementsDeg

function Orbit(args...; kwargs...)
    @warn "Orbit is deprecated in favour of KeplerianElements"
    return KeplerianElements(args...; kwrags...)
end
export Orbit

# Better printing
Base.show(io::IO, ::MIME"text/plain", elem::KeplerianElements) = print(
    io, """
        $(typeof(elem))
        ─────────────────────────
        a   [au ] = $(round(elem.a,sigdigits=3)) 
        i   [°  ] = $(round(rad2deg(elem.i),sigdigits=3))
        e         = $(round(elem.e,sigdigits=3))
        τ         = $(round(elem.τ,sigdigits=3))
        μ   [M⊙ ] = $(round(elem.μ,sigdigits=3)) 
        ω   [°  ] = $(round(rad2deg(elem.ω),sigdigits=3))
        Ω   [°  ] = $(round(rad2deg(elem.Ω),sigdigits=3))
        plx [mas] = $(round(elem.plx,sigdigits=3)) 
        ──────────────────────────
        period      [yrs ] : $(round(period(elem)/year2days,digits=1)) 
        distance    [pc  ] : $(round(distance(elem),digits=1)) 
        mean motion [°/yr] : $(round(rad2deg(meanmotion(elem)),sigdigits=3)) 
        ──────────────────────────
        """)
Base.show(io::IO, elem::KeplerianElements) = print(io,
    "KeplerianElements($(round(elem.a,sigdigits=3)), $(round(elem.i,sigdigits=3)), $(round(elem.e,sigdigits=3)), "*
    "$(round(elem.τ,sigdigits=3)), $(round(elem.μ,sigdigits=3)), $(round(elem.ω,sigdigits=3)), "*
    "$(round(elem.Ω,sigdigits=3)), $(round(elem.plx,sigdigits=3)))"
)

Base.show(io::IO, ::MIME"text/html", elem::KeplerianElements) = print(
    io, """
        <table style="font-family:monospace; text-align: right">
        <tr><th colspan=3 style="font-family:sans-serif; text-align: left">$(typeof(elem))</th></tr>
        <tr><td rowspan=8>Input</td><td>a   [au] =</td> <td>$(round(elem.a,sigdigits=3))</td></tr>
        <tr><td>i   [°] = </td><td>$(round(rad2deg(elem.i),sigdigits=3))</td></tr>
        <tr><td>e         = </td><td>$(round(elem.e,sigdigits=3))</td></tr>
        <tr><td>τ         = </td><td>$(round(elem.τ,sigdigits=3))</td></tr>
        <tr><td>μ   [M⊙] = </td><td>$(round(elem.μ,sigdigits=3)) </td></tr>
        <tr><td>ω   [°] = </td><td>$(round(rad2deg(elem.ω),sigdigits=3))</td></tr>
        <tr><td>Ω   [°] = </td><td>$(round(rad2deg(elem.Ω),sigdigits=3))</td></tr>
        <tr><td>plx [mas] = </td><td>$(round(elem.plx,sigdigits=3)) </td></tr>
        <tr><td rowspan=3>Computed</td><td>period      [yrs] : </td><td>$(round(period(elem)/DirectOrbits.year2days,digits=1)) </td></tr>
        <tr><td>distance    [pc] : </td><td>$(round(distance(elem),digits=1)) </td></tr>
        <tr><td>mean motion [°/yr] : </td><td>$(round(rad2deg(DirectOrbits.meanmotion(elem)),sigdigits=3)) </td></tr>
        </table>
        """)

import Base: length, iterate
length(::AbstractElements) = 1
iterate(elem::AbstractElements) = (elem, nothing)
iterate(::AbstractElements, ::Nothing) = nothing


"""
    period(elem)

Period of an orbit in days.
"""
period(elem::KeplerianElements) = elem.T
export period

"""
    distance(elem)

Distance to the system in parsecs.
"""
distance(elem::KeplerianElements) = elem.dist/pc2au
export distance

"""
    meanmotion(elem)

Mean motion, radians per year.
"""
meanmotion(elem::KeplerianElements) = elem.n

# Data type template for results of kep2cart
const template = ComponentVector(x=1.,y=1.,z=1.,ẋ=1.,ẏ=1.,ż=1.)
const template_axes = ComponentArrays.getaxes(template)
const Position = typeof(template)

"""
    kep2cart(elements, t)

Given an set of orbital elements with a time `t` in days, get the position and
velocity of the secondary body (e.g. planet around a star).

This will output a static component array with the following properties:
 - `x`: δ right ascension in milliarcseconds
 - `y`: δ declination in milliarcseconds
 - `z`: δ distance along the line of sight in "milliarcseconds" for consistency
 - `ẋ`: right ascension proper motion anomaly in milliarcseconds/year
 - `ẏ`: declination proper motion anomaly in milliarcseconds/year
 - `ż`: radial velocity of the *componanion* in meters/second.

You can access the properties either by name `.x` or by index `[1]`. There are helper
functions to calculate each of these properties individually, but if you need more than
one it is most efficient to calculate them in one go.

`radvel` can optionally accept the masses of the central and secondary
bodies to calculate the impact of the secondary body on radial velocity of the star,
instead of the radial velocity of the secondary body itself.

Note: these calculations use the small angle approximation, so are only accurate when 
the star is much further way from the observer than the companion is from the star.

See also: `kep2cart_ν`, `projectedseparation`, `raoff`, `decoff`, `radvel`, `propmotionanom`.
"""
@inline function kep2cart(elem::KeplerianElements{T}, t; tref=58849) where T
    T2 = promote_type(T, typeof(t))
    
    # Compute mean anomaly
    tₚ = elem.τ*period(elem) + tref
    
    MA = meanmotion(elem)/convert(T2, year2days) * (t - tₚ)

    if !isfinite(MA)
        MA = zero(typeof(MA))
        @warn "non-finite mean anomaly" maxlog=50
    end 

    # Compute eccentric anomaly
    # This uses the kepler_solver function from AstroLib.
    # It's by far the fastest function for solving Kepler's
    # equation that I have tested.
    EA = _kepler_solver_inline(MA, elem.e)

   
    if !isfinite(EA)
        EA = MA
        @warn "non-finite eccentric anomaly" elem.e maxlog=50
    end
    
    # Calculate true anomaly
    ν = convert(T2,2)*atan(elem.ν_fact*tan(EA/convert(T2,2)))
    # ν = convert(T2,2)*elem.ν_fact*(EA/convert(T2,2))

    if !isfinite(ν)
        ν = zero(typeof(ν))
        @warn "non-finite true anomaly" maxlog=50
    end
    

    # New radius.
    # This is the semi-major axis, modified by the eccentricity. Units of AU.
    r = elem.a*(one(T2)-elem.e*cos(EA))

    
    if !isfinite(r)
        r = zero(typeof(r))
        @warn "non-finite radius"
    end
    
    
    # Project back into Cartesian coordinates (AU).
    sin_ω_ν, cos_ω_ν = sincos(elem.ω+ν)
    xₐᵤ = r*(elem.sin_Ω*cos_ω_ν + elem.cos_Ω*sin_ω_ν*elem.cos_i)
    yₐᵤ = r*(elem.cos_Ω*cos_ω_ν - elem.sin_Ω*sin_ω_ν*elem.cos_i)
    zₐᵤ = r*(elem.sin_i*sin_ω_ν)

    # Radial velocity
    p = elem.a*(1-elem.e^2) 
    h = sqrt(elem.μ*p) # Specific angular momentum
    
    # Note: use the small angle approximation since arctangent is relatively slow.
    dist⁻¹ = 1/elem.dist
    xᵣ = xₐᵤ*dist⁻¹ # atan(xₐᵤ, elem.dist)
    yᵣ = yₐᵤ*dist⁻¹ # atan(yₐᵤ, elem.dist)
    zᵣ = zₐᵤ*dist⁻¹ # atan(zₐᵤ, elem.dist)

    xₘₐₛ = xᵣ * rad2as*oftype(xᵣ,1e3)
    yₘₐₛ = yᵣ * rad2as*oftype(yᵣ,1e3)
    zₘₐₛ = zᵣ * rad2as*oftype(zᵣ,1e3)

    # Factor out common sub-expressions
    r⁻¹ = 1/r
    A = h * elem.e / p * r⁻¹ * sin(ν)
    h_r = h * r⁻¹
    cos_ω_ν_cos_i = cos_ω_ν*elem.cos_i

    # TODO: figure out units from first principles.
    ẋₐᵤ = xₐᵤ * A - h_r*(elem.cos_Ω*sin_ω_ν + elem.sin_Ω*cos_ω_ν_cos_i)
    ẏₐᵤ = yₐᵤ * A - h_r*(elem.sin_Ω*sin_ω_ν - elem.cos_Ω*cos_ω_ν_cos_i)
    żₐᵤ = zₐᵤ * A + h_r*elem.sin_i*cos_ω_ν
    
    # We want radial velocity in m/s, and the tangential velocities
    # in mas/year. The
    # ẋₖₘₛ = ẋₐᵤ * 29780
    # ẏₖₘₛ = ẏₐᵤ * 29780
    żₖₘₛ = żₐᵤ * 29780

    # Note: use the small angle approximation since arctangent is relatively slow.
    ẋᵣ = ẋₐᵤ*dist⁻¹ # atan(ẋₐᵤ, elem.dist)
    ẏᵣ = ẏₐᵤ*dist⁻¹ # atan(ẏₐᵤ, elem.dist)
    # żᵣ = żₐᵤ*dist⁻¹ # atan(żₐᵤ, elem.dist)

    ẋₘₐₛₐ = ẋᵣ * rad2as*oftype(ẋᵣ,1e3) * 2π
    ẏₘₐₛₐ = ẏᵣ * rad2as*oftype(ẏᵣ,1e3) * 2π
    # zₘₐₛ = zᵣ * rad2as*oftype(zᵣ,1e3)

    return ComponentVector(SVector(xₘₐₛ, yₘₐₛ, zₘₐₛ, ẋₘₐₛₐ, ẏₘₐₛₐ, żₖₘₛ), template_axes)
    # return (;
    #     x=xₘₐₛ,
    #     y=yₘₐₛ,
    #     z=zₘₐₛ,
    #     ẋ=ẋₘₐₛₐ,
    #     ẏ=ẏₘₐₛₐ,
    #     ż=żₖₘₛ
    # )
end
export kep2cart

# Kep2cart, only it directly accepts the true anomaly.
# A key use for this is tracing out an orbit for a plot, where
# you want roughly equal spacing of points by angle, rather
# than in time (bunching up at apoapsis and not enough at periapsis)
function kep2cart_ν(elem::KeplerianElements, ν)
    
    # Semi-latus of rectum    
    p = elem.a*(1-elem.e^2) 
    r = p/(1+elem.e*cos(ν))
    
    # Project back into Cartesian coordinates (AU).
    sin_ω_ν, cos_ω_ν = sincos(elem.ω+ν)
    xₐᵤ = r*(elem.sin_Ω*cos_ω_ν + elem.cos_Ω*sin_ω_ν*elem.cos_i)
    yₐᵤ = r*(elem.cos_Ω*cos_ω_ν - elem.sin_Ω*sin_ω_ν*elem.cos_i)
    zₐᵤ = r*(elem.sin_i*sin_ω_ν)

    # Radial velocity
    h = sqrt(elem.μ*p) # Specific angular momentum
    
    dist⁻¹ = 1/elem.dist
    # Note: use the small angle approximation since arctangent is relatively slow.
    xᵣ = xₐᵤ*dist⁻¹ # atan(xₐᵤ, elem.dist)
    yᵣ = yₐᵤ*dist⁻¹ # atan(yₐᵤ, elem.dist)
    zᵣ = zₐᵤ*dist⁻¹ # atan(zₐᵤ, elem.dist)

    xₘₐₛ = xᵣ * rad2as*oftype(xᵣ,1e3)
    yₘₐₛ = yᵣ * rad2as*oftype(yᵣ,1e3)
    zₘₐₛ = zᵣ * rad2as*oftype(zᵣ,1e3)

    # Factor out common sub-expressions
    A = h * elem.e / (r*p) * sin(ν)
    h_r = h / r

    # TODO: figure out units from first principles.
    ẋₐᵤ = xₐᵤ * A - h_r*(elem.cos_Ω*sin_ω_ν + elem.sin_Ω*cos_ω_ν*elem.cos_i)
    ẏₐᵤ = yₐᵤ * A - h_r*(elem.sin_Ω*sin_ω_ν - elem.cos_Ω*cos_ω_ν*elem.cos_i)
    żₐᵤ = zₐᵤ * A + h_r*elem.sin_i*cos_ω_ν
    
    # We want radial velocity in m/s, and the tangential velocities
    # in mas/year. The
    # ẋₖₘₛ = ẋₐᵤ * 29780
    # ẏₖₘₛ = ẏₐᵤ * 29780
    żₖₘₛ = żₐᵤ * 29780

    # Note: use the small angle approximation since arctangent is relatively slow.
    ẋᵣ = ẋₐᵤ*dist⁻¹ # atan(ẋₐᵤ, elem.dist)
    ẏᵣ = ẏₐᵤ*dist⁻¹ # atan(ẏₐᵤ, elem.dist)
    # żᵣ = żₐᵤ*dist⁻¹ # atan(żₐᵤ, elem.dist)

    # TODO: investigate source of 2pi factor
    ẋₘₐₛₐ = ẋᵣ * rad2as*oftype(ẋᵣ,1e3) * 2π
    ẏₘₐₛₐ = ẏᵣ * rad2as*oftype(ẏᵣ,1e3) * 2π
    # zₘₐₛ = zᵣ * rad2as*oftype(zᵣ,1e3)

    return ComponentVector(SVector(xₘₐₛ, yₘₐₛ, zₘₐₛ, ẋₘₐₛₐ, ẏₘₐₛₐ, żₖₘₛ), template_axes)
end


"""
    raoff(elements, t)

Get the offset from the central body in Right Ascention in
milliarcseconds at some time `t` in days.
"""
function raoff(elements::AbstractElements, t)
    return kep2cart(elements, t).x
end
export raoff

"""
    decoff(elements, t)

Get the offset from the central body in Declination in
milliarcseconds at some time `t` in days.
"""
function decoff(elements::AbstractElements, t)
    return kep2cart(elements, t).y
end
export decoff

"""
    losoff(elements, t)

Get the offset from the central body in the line of sight towards
the system at time `t` in days, also in milliarcseconds. Of course, we can't observe this
as an angle, but we use the same units for consistency.
"""
function losoff(elements::AbstractElements, t)
    return kep2cart(elements, t).z
end
export losoff

"""
    radvel(elements, t)

Get the radial velocity of the *planet* along the line of sight
at the time `t` in days, in units of m/s.
"""
function radvel(elements::AbstractElements, t)
    return kep2cart(elements, t).ż
end
export radvel

"""
    radvel(elements, t, M_star, M_planet)

Get the radial velocity of the *star* along the line of sight
at the time `t` in days, in units of m/s.
The mass of the star and planet must have consistent units.
"""
function radvel(elements::AbstractElements, t, M_planet)
    M_star = elements.μ
    v_planet =  kep2cart(elements, t).ż
    v_star = v_planet * M_planet / M_star
    return v_star
end
export radvel

"""
    propmotionanom(elements, t, M_star, M_planet)

Calculate the instantenous proper motion anomaly on a star due to
an orbiting companion.
"""
function propmotionanom(elements::AbstractElements, t, M_planet)
    M_star = elements.μ
    o = kep2cart(elements, t)
    Δμ_planet = -SVector(o.ẋ, o.ẏ) # milliarcseconds per year
    Δμ_star = Δμ_planet * M_planet / (M_star + M_planet)
    return Δμ_star
end
function propmotionanom(elements::AbstractElements, t)
    o = kep2cart(elements, t)
    Δμ_planet = -SVector(o.ẋ, o.ẏ) # milliarcseconds per year
    return Δμ_planet
end
export propmotionanom

"""
    posangle(elements, t)

Calculate the instantenous proper motion anomaly on a star due to
an orbiting companion.
"""
function posangle(elements::AbstractElements, t, M_star, M_planet)
    o = kep2cart(elements, t)
    return atan(o.y,o.x)
end
export propmotionanom

"""
    projectedseparation(elements, t)

Projected separation in mas from the central body at time t (days).
"""
function projectedseparation(elements::AbstractElements, t)
    x,y,z = kep2cart(elements,t)
    return sqrt(x^2 + y^2 + z^2)
end
export projectedseparation

using RecipesBase
@recipe function f(elem::AbstractElements)
    νs = range(-π, π, length=90)
    coords = kep2cart_ν.(elem, νs)
    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]

    # We almost always want to see spatial coordinates with equal step sizes
    aspect_ratio --> 1
    # And we almost always want to reverse the RA coordinate to match how we
    # see it in the sky.
    xflip --> true

    return xs, ys
end

@recipe function f(elems::AbstractArray{<:AbstractElements})

    # Step through true anomaly instead of time.
    # This produces far nicer looking plots, especially if
    # the orbits in question vary significantly in period
    # or are eccentric
    νs = range(-π, π, length=90)
    coords = kep2cart_ν.(elems, νs')

    xs = [c[1] for c in coords]'
    ys = [c[2] for c in coords]'

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

    seriesalpha --> 30/length(elems)


    return xs, ys
end

include("Transformation.jl")
include("Time.jl")

# The following function is taken directly from AstroLib.jl
# We  remove one invariant check we handle elsewhere and also
# force inlining for about a 5% speedup.
# We also supply analytic gradients for use in autodiff packages.
@inline function _kepler_solver_inline(_M::Real, e::Real)
    # We already handle this invariant
    # @assert 0 <= e <= 1 "eccentricity must be in the range [0, 1]"
    # M must be in the range [-pi, pi], see Markley (1995), page 2.
    # M = rem2pi(_M, RoundNearest)
    M = rem2pi_safe(_M)
    T = float(promote_type(typeof(M), typeof(e)))
    if iszero(M) || iszero(e)
        return T(M)
    end
    pi2 = abs2(T(pi))
    # equation (20)
    α = (3 * pi2 + 8 * (pi2 - pi * abs(M)) / (5 * (1 + e)))/(pi2 - 6)
    # equation (5)
    d = 3 * (1 - e) + α * e
    # equation (9)
    q = 2 * α * d * (1 - e) - M * M
    # equation (10)
    r = 3 * α * d * (d - 1 + e) * M + M * M * M
    # equation (14)
    w = cbrt(abs2(abs(r) + sqrt(q * q * q + r * r)))
    # equation (15)
    E1 = (2 * r * w / @evalpoly(w, q * q, q, 1) + M)/d
    # equation (26) & equation (27)
    f2, f3 = e .* sincos(E1)
    # equation (21)
    f0 = E1 - f2 - M
    # equation (25)
    f1 = 1 - f3
    # equation (22)
    δ3 = -f0 / (f1 - f0 * f2 / (2 * f1))
    # equation (23)
    δ4 = -f0 / @evalpoly(δ3, f1, f2 / 2, f3 / 6)
    # equations (24) and (28)
    δ5 = -f0 / @evalpoly(δ4, f1, f2 / 2, f3 / 6, - f2 / 24)
    return E1 + δ5 # equation 29
end

# Using implicit differentiation, I found that the derivatives of eccentric anomaly
# have closed form solutions once the primal value is known. 
# By providing those here, upstream automatic differentiation libraries will be able
# to efficiently diff through Kepler's equation.
# using ChainRulesCore
# @scalar_rule _kepler_solver_inline(MA, e) @setup(u = 1 - e*cos(Ω)) (1 / u,sin(Ω) / u)

# We try to support symbolic manipulation using Symbolics.jl, but it's
# not reasonable to use `remp2pi` on a symbolic variable.
# We therefore have a special fallback method for that case. We 
# define it when both packages get loaded by the user using Requires.
@inline rem2pi_safe(x) = rem2pi(x, RoundNearest)

using ChainRulesCore
# @scalar_rule _kepler_solver_inline(MA, e) @setup(u = 1 - e*cos(Ω)) (1 / u,sin(Ω) / u)

# Define a scale rule to allow Zygote to diff through rem2pi
@scalar_rule rem2pi_safe(x) x


using Requires
function __init__()
    @require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" begin
        @inline rem2pi_safe(x::Symbolics.Num) = x
    end
end

end # module