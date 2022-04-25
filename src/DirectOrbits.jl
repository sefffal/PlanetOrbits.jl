"""
# DirectOrbits
A package for calculating orbits in the context of direct imaging,
astrometry, and radial velocity.
"""

module DirectOrbits

# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

using LinearAlgebra
using StaticArrays

# ----------------------------------------------------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------------------------------------------------

# radians <-> milliarcseconds
const rad2mas = 2.06264806e8
const mas2rad = 4.848136816903219e-9

# radians <-> arcseconds 
const rad2as = 206265
const as2rad = 4.848132257047972e-6

# parsecs <-> astronomical units
const pc2au = 206265
const au2pc = 4.848132257047972e-6

# astronomical units <-> metres
const au2m = 1.495978707e11
const m2au = 6.684587122268445e-12

# years <-> days
const year2day = 365.2422
const day2year = 2.737909255830788e-3

# years <-> seconds
const year2sec = 3.1556926e7
const sec2year = 3.168876461541279e-8

# days <-> seconds
const day2sec = 86400
const sec2day = 1.1574074074074073e-5

# ----------------------------------------------------------------------------------------------------------------------
# KeplerianElements
# ----------------------------------------------------------------------------------------------------------------------

abstract type AbstractOrbit end
export AbstractOrbit

"""
    KeplerianElements(
        a, # semi-major axis [AU]
        e, # eccentricity
        i, # inclination [rad]
        ω, # argument of periapsis [rad]
        Ω, # longitude of ascending node [rad]
        τ, # epoch of periastron passage at MJD=0
        M, # mass of primary [M⊙]
        plx, # parallax [mas]; defines the distance to the primary
    )

Represents the Keplerian elements of a secondary body orbiting a primary.
Values can be specified by keyword argument or named tuple for convenience.

See also `KeplerianElementsDeg` for a convenience constructor accepting
units of degrees instead of radians for `i`, `ω`, and `Ω`.
"""
struct KeplerianElements{T<:Number} <: AbstractOrbit

    # Orbital properties
    a::T
    e::T
    i::T
    ω::T
    Ω::T
    τ::T
    M::T
    plx::T

    # Physical constants
    dist::T
    T::T
    n::T
    ν_fact::T
    p::T

    # Geometric factors
    cosi::T
    sini::T
    cosω::T
    sinω::T
    cosΩ::T
    sinΩ::T
    ecosω::T
    esinω::T
    cosi_cosΩ::T
    cosi_sinΩ::T

    # Semiamplitudes
    J::T
    K::T
    A::T

    # Inner constructor to enforce invariants and pre-calculate
    # constants from the orbital elements
    function KeplerianElements(a, e, i, ω, Ω, τ, M, plx)

        # Enforce invariants on user parameters
        a = max(a, zero(a))
        e = max(zero(e), min(e, one(e)))
        τ = mod(τ, one(τ))
        M = max(M, zero(M))
        plx = max(plx, zero(plx))

        # Pre-calculate factors to be re-used by orbitsolve
        # Physical constants of system and orbit
        dist = 1000/plx * pc2au # distance [AU]
        rootacubeoverm = √(a^3/M)
        periodyrs = rootacubeoverm
        period = periodyrs * year2day # period [days]
        n = 2π/√(a^3/M) # mean motion
        ν_fact = √((1 + e)/(1 - e)) # true anomaly prefactor
        oneminusesq = (1 - e^2)
        p = a*oneminusesq # semi-latus rectum [AU]

        # Get type of parameters
        T = promote_type(
            typeof(a), typeof(e), typeof(i), typeof(ω),
            typeof(Ω), typeof(τ), typeof(M), typeof(plx),
        )

        # The user might pass in integers, but it makes no sense to do these
        # calculations on integers. Assume they mean to use floats.
        if T <: Integer
            T = promote_type(T, Float64)
        end

        # Geometric factors involving rotation angles
        sini, cosi = sincos(i)
        sinω, cosω = sincos(ω)
        sinΩ, cosΩ = sincos(Ω)
        ecosω = e*cosω
        esinω = e*sinω
        cosi_cosΩ = cosi*cosΩ
        cosi_sinΩ = cosi*sinΩ

        # Velocity and acceleration semiamplitudes
        J = ((2π*a)/periodyrs) / √oneminusesq # horizontal velocity semiamplitude [AU/year]
        K = J*au2m*sec2year*sini # radial velocity semiamplitude [m/s]
        A = ((4π^2 * a)/periodyrs^2) / oneminusesq^2 # horizontal acceleration semiamplitude [AU/year^2]

        new{T}(
            # Passed parameters that define the elements
            a, e, i, ω, Ω, τ, M, plx,
            # Cached calcuations
            dist, period, n, ν_fact, p,
            # Geometric factors
            cosi, sini, cosω, sinω, cosΩ, sinΩ, ecosω, esinω, cosi_cosΩ, cosi_sinΩ,
            # Semiamplitudes
            J, K, A
        )
    end
end

# Allow arguments to be specified by keyword
KeplerianElements(;a, e, i, ω, Ω, τ, M, plx) = KeplerianElements(a, e, i, ω, Ω, τ, M, plx)
# Allow arguments to be specified by named tuple
KeplerianElements(nt) = KeplerianElements(nt.a, nt.e, nt.i, nt.ω, nt.Ω, nt.τ, nt.M, nt.plx)
export KeplerianElements

"""
    KeplerianElementsDeg(a, e, i, ω, Ω, τ, M, plx)

A convenience function for constructing KeplerianElements where
`i`, `ω`, and `Ω` are provided in units of degrees instead of radians.
"""
KeplerianElementsDeg(a, e, i, ω, Ω, τ, M, plx) = KeplerianElements(a, e, deg2rad(i), deg2rad(ω), deg2rad(Ω), τ, M, plx)
KeplerianElementsDeg(;a, e, i, ω, Ω, τ, M, plx) = KeplerianElementsDeg(a, e, i, ω, Ω, τ, M, plx)
KeplerianElementsDeg(nt) = KeplerianElementsDeg(nt.a, nt.e, nt.i, nt.ω, nt.Ω, nt.τ, nt.M, nt.plx)
export KeplerianElementsDeg

"""
    astuple(elements)

Return the parameters of a KeplerianElements value as a tuple.
"""
function astuple(elem::KeplerianElements)
    return (;elem.a, elem.e, elem.i, elem.ω, elem.Ω, elem.τ, elem.M, elem.plx)
end
export astuple

# Pretty printing
Base.show(io::IO, ::MIME"text/plain", elem::KeplerianElements) = print(
    io, """
        $(typeof(elem))
        ─────────────────────────
        a   [au ] = $(round(elem.a, sigdigits=3))
        e         = $(round(elem.e, sigdigits=3))
        i   [°  ] = $(round(rad2deg(elem.i), sigdigits=3))
        ω   [°  ] = $(round(rad2deg(elem.ω), sigdigits=3))
        Ω   [°  ] = $(round(rad2deg(elem.Ω), sigdigits=3))
        τ         = $(round(elem.τ, sigdigits=3))
        M   [M⊙ ] = $(round(elem.M, sigdigits=3)) 
        plx [mas] = $(round(elem.plx, sigdigits=3)) 
        ──────────────────────────
        period      [yrs ] : $(round(period(elem)*day2year, digits=1)) 
        distance    [pc  ] : $(round(distance(elem), digits=1)) 
        mean motion [°/yr] : $(round(rad2deg(meanmotion(elem)), sigdigits=3)) 
        ──────────────────────────
        """
)

Base.show(io::IO, elem::KeplerianElements) = print(io,
    "KeplerianElements($(round(elem.a, sigdigits=3)), $(round(elem.e, sigdigits=3)), $(round(elem.i, sigdigits=3)), "*
    "$(round(elem.ω, sigdigits=3)), $(round(elem.Ω, sigdigits=3)), $(round(elem.τ, sigdigits=3)), "*
    "$(round(elem.M, sigdigits=3)), $(round(elem.plx, sigdigits=3)))"
)

# Pretty printing in notebooks as HTML
Base.show(io::IO, ::MIME"text/html", elem::KeplerianElements) = print(
    io, """
        <table style="font-family:monospace; text-align: right">
        <tr><th colspan=3 style="font-family:sans-serif; text-align: left">$(typeof(elem))</th></tr>
        <tr><td rowspan=8>Input</td><td>a   [au] =</td> <td>$(round(elem.a, sigdigits=3))</td></tr>
        <tr><td>e         = </td><td>$(round(elem.e, sigdigits=3))</td></tr>
        <tr><td>i   [°] = </td><td>$(round(rad2deg(elem.i), sigdigits=3))</td></tr>
        <tr><td>ω   [°] = </td><td>$(round(rad2deg(elem.ω), sigdigits=3))</td></tr>
        <tr><td>Ω   [°] = </td><td>$(round(rad2deg(elem.Ω), sigdigits=3))</td></tr>
        <tr><td>τ         = </td><td>$(round(elem.τ, sigdigits=3))</td></tr>
        <tr><td>M   [M⊙] = </td><td>$(round(elem.M, sigdigits=3)) </td></tr>
        <tr><td>plx [mas] = </td><td>$(round(elem.plx, sigdigits=3)) </td></tr>
        <tr><td rowspan=3>Computed</td><td>period      [yrs] : </td><td>$(round(period(elem)*DirectOrbits.day2year, digits=1)) </td></tr>
        <tr><td>distance    [pc] : </td><td>$(round(distance(elem), digits=1)) </td></tr>
        <tr><td>mean motion [°/yr] : </td><td>$(round(rad2deg(DirectOrbits.meanmotion(elem)), sigdigits=3)) </td></tr>
        </table>
        """
)

# Define iterate and length = 1 so that we can broadcast over elements.
Base.length(::AbstractOrbit) = 1
Base.iterate(elem::AbstractOrbit) = (elem, nothing)
Base.iterate(::AbstractOrbit, ::Nothing) = nothing

# ----------------------------------------------------------------------------------------------------------------------
# Orbit Solutions
# ----------------------------------------------------------------------------------------------------------------------
abstract type AbstractOrbitSolution end

"""
    OrbitSolutionKeplerian(
        x, # δ right ascension [mas]
        y, # δ declination [mas]
        ẋ, # right ascension proper motion anomaly [mas/year]
        ẏ, # declination proper motion anomaly [mas/year]
        ż, # radial velocity of the *secondary* [m/s]
        ẍ, # right ascension acceleration [mas/year^2]
        ÿ, # declination acceleration [mas/year^2]
        elem, # KeplerianElements representing the orbit of this body 
    )

Represents the secondary's position on the sky in terms of offset from
the primary, its velocity and acceleration on the sky, and its radial velocity.
Conceptually, this is a KeplerianElements evaluated to some position.
"""
struct OrbitSolutionKeplerian{T<:Number,TEl<:KeplerianElements} <: AbstractOrbitSolution
    elem::TEl
    ν::T
    sinν_ω::T
    cosν_ω::T
    ecosν::T
    r::T
    function OrbitSolutionKeplerian(elem, ν, sinν_ω, cosν_ω, ecosν, r)
        promoted = promote(ν, sinν_ω, cosν_ω, ecosν, r)
        return new{eltype(promoted),typeof(elem)}(elem, promoted...)
    end
end
export OrbitSolutionKeplerian

# Printing
# Base.show(io::IO, os::AbstractOrbitSolution) = print(io,
#     "AbstractOrbitSolution(x = $(round(os.x, sigdigits=3)), y = $(round(os.y, sigdigits=3)), "*
#     "ẋ = $(round(os.ẋ, sigdigits=3)), ẏ = $(round(os.ẏ, sigdigits=3)), ż = $(round(os.ż, sigdigits=3)), "*
#     "ẍ = $(round(os.ẍ, sigdigits=3)), ÿ = $(round(os.ÿ, sigdigits=3)))"
# )

# ----------------------------------------------------------------------------------------------------------------------
# System Properties
# ----------------------------------------------------------------------------------------------------------------------

"""
    period(elem)

Period of an orbit [days].
"""
period(elem::KeplerianElements) = elem.T
export period

"""
    distance(elem)

Distance to the system [pc].
"""
distance(elem::KeplerianElements) = elem.dist*au2pc
export distance

"""
    meanmotion(elem)

Mean motion [rad/year].
"""
meanmotion(elem::KeplerianElements) = elem.n
export meanmotion 

"""
   periastron(elements, tref=58849)

Compute the MJD of periastron passage most recently after the reference epoch tref.
N.B. mjd of 58849 = 2020-01-01
"""
function periastron(elem::AbstractOrbit, tref=58849)
    tₚ = elem.τ*period(elem) + tref
    return tₚ
end
export periastron

"""
    semiamplitude(elem)

Radial velocity semiamplitude [m/s].
"""
semiamplitude(elem::KeplerianElements) = elem.K
export semiamplitude 

# ----------------------------------------------------------------------------------------------------------------------
# Solve Orbit in Cartesian Coordinates
# ----------------------------------------------------------------------------------------------------------------------

"""
    orbitsolve_ν(elem, ν)

Solve a keplerian orbit from a given true anomaly [rad].
See orbitsolve for the same function accepting a given time.
"""
function orbitsolve_ν(elem::KeplerianElements, ν)
    sinν_ω, cosν_ω = sincos(elem.ω + ν)
    ecosν = elem.e*cos(ν)
    r = elem.p/(1 + ecosν)
    return OrbitSolutionKeplerian(elem, ν, sinν_ω, cosν_ω, ecosν, r)
end

"""
    orbitsolve(elements, t)

Given a set of orbital elements with a time `t` in days, get the position and
velocity of the secondary body (e.g. planet around a star).

This will output an `AbstractOrbitSolution` struct with the following properties:
 - `x`: δ right ascension [mas]
 - `y`: δ declination [mas]
 - `ẋ`: right ascension proper motion anomaly [mas/year]
 - `ẏ`: declination proper motion anomaly [mas/year]
 - `ż`: radial velocity of the *secondary* [m/s]
 - `ẍ`: right ascension acceleration [mas/year^2]
 - `ÿ`: declination acceleration [mas/year^2]

You can access the properties by name `.x`. There are helper functions to
calculate each of these properties individually, but if you need more than
one it is most efficient to calculate them in one go.

`radvel` can optionally accept the mass of the primary to calculate the impact
of the secondary body on radial velocity of the primary, instead of the radial
velocity of the secondary body itself.

Note: these calculations use the small angle approximation, so are only accurate when 
the star is much further way from the observer than the secondary is from the primary.

See also: `orbitsolve_ν`, `projectedseparation`, `raoff`, `decoff`, `radvel`, `propmotionanom`.
"""
function orbitsolve(elem::KeplerianElements{T}, t; tref=58849) where T
    
    # Epoch of periastron passage
    tₚ = periastron(elem, tref)

    # Mean anomaly    
    MA = meanmotion(elem)/oftype(t, year2day) * (t - tₚ)

    # Compute eccentric anomaly
    EA = kepler_solver(MA, elem.e)
    
    # Calculate true anomaly
    ν = 2*atan(elem.ν_fact*tan(EA/2))

    return orbitsolve_ν(elem, ν)
end
export orbitsolve

# ----------------------------------------------------------------------------------------------------------------------
# Orbital Position and Motion
# ----------------------------------------------------------------------------------------------------------------------

"""
    raoff(elem, t)

Get the offset [mas] from the primary body in Right Ascension
at the time `t` [days].

    raoff(o)

Get the offset [mas] from the primary body in Right Ascension 
from an instance of `AbstractOrbitSolution`.
"""
function raoff(o::AbstractOrbitSolution)
    xcart = o.r*(o.cosν_ω*o.elem.sinΩ + o.sinν_ω*o.elem.cosi*o.elem.cosΩ) # [AU]
    cart2angle = rad2as*oftype(xcart, 1e3)/o.elem.dist
    xang = xcart*cart2angle # [mas]
    return xang
end
export raoff

"""
    decoff(elem, t)

Get the offset [mas] from the primary body in Declination
at the time `t` [days].

    decoff(elem, t)

Get the offset [mas] from the primary body in Declination
from an instance of `AbstractOrbitSolution`.
"""
function decoff(o::AbstractOrbitSolution)
    ycart = o.r*(o.cosν_ω*o.elem.cosΩ - o.sinν_ω*o.elem.cosi*o.elem.sinΩ) # [AU]
    cart2angle = rad2as*oftype(ycart, 1e3)/o.elem.dist
    yang = ycart*cart2angle # [mas]
    return yang
end
export decoff

"""
    posangle(elem, t)

Calculate the position angle [rad] of the secondary about its primary
from our perspective at the time `t` [days].

    posangle(o)

Calculate the position angle [rad] of the secondary about its primary
from our perspective from an instance of `AbstractOrbitSolution`.
"""
function posangle(o::AbstractOrbitSolution)
    x = raoff(o)
    y = decoff(o)
    return atan(x, y) # Note: the order of these arguments is *correct* in our conventions
end
export posangle

"""
    projectedseparation(elem, t)

Calculate the projected separation [mas] of the secondary from its
primary at the time `t` [days].

    projectedseparation(o)

Calculate the projected separation [mas] of the secondary from its
primary from an instance of `AbstractOrbitSolution`.
"""
function projectedseparation(o::AbstractOrbitSolution)
    x = raoff(o)
    y = decoff(o)
    return sqrt(x^2 + y^2)
end
export projectedseparation

"""
    propmotionanom(elem, t)

Get the instantaneous proper motion anomaly [mas/year] of
the *secondary* at the time `t` [days].

    propmotionanom(o)

Get the instantaneous proper motion anomaly [mas/year] of
the *secondary* from an instance of `AbstractOrbitSolution`.
"""
function propmotionanom(o::AbstractOrbitSolution)
    Δμ_planet = SVector(pmra(o), pmdec(o))
    return Δμ_planet
end

function pmra(o::AbstractOrbitSolution)
    ẋcart = o.elem.J*(o.elem.cosi_cosΩ*(o.cosν_ω + o.elem.ecosω) - o.elem.sinΩ*(o.sinν_ω + o.elem.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẋcart, 1e3)/o.elem.dist
    ẋang = ẋcart*cart2angle # [mas/year]
    return ẋang
end

function pmdec(o::AbstractOrbitSolution)
    ẏcart = -o.elem.J*(o.elem.cosi_sinΩ*(o.cosν_ω + o.elem.ecosω) + o.elem.cosΩ*(o.sinν_ω + o.elem.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẏcart, 1e3)/o.elem.dist
    ẏang = ẏcart*cart2angle # [mas/year]
    return ẏang
end

export pmra, pmdec

"""
    propmotionanom(elem, t, M_planet)

Get the instantaneous proper motion anomaly [mas/year] of 
the *primary* in at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.
"""
function propmotionanom(o::AbstractOrbitSolution, M_planet)
    M_star = o.elem.M
    Δμ_planet = propmotionanom(o)
    Δμ_star = -(M_planet/(M_star + M_planet))*Δμ_planet
    return Δμ_star
end
export propmotionanom

"""
    radvel(elem, t)

Get the radial velocity [m/s] of the *secondary* along the
line of sight at the time `t` [days].

    radvel(o)

Get the radial velocity [m/s] of the *secondary* along the
line of sight from an instance of `AbstractOrbitSolution`.
"""
function radvel(o::AbstractOrbitSolution)
    żcart = o.elem.K*(o.cosν_ω + o.elem.ecosω) # [m/s]
    return żcart
end

"""
    radvel(elem, t, M_planet)

Get the radial velocity [m/s] of the *primary* along the
line of sight at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.

radvel(elem, t, M_planet)

Get the radial velocity [m/s] of the *primary* along the
line of sight from an `AbstractOrbitSolution`. The units of `M_planet`
and `elem.M` must match.
"""
function radvel(o::AbstractOrbitSolution, M_planet)
    M_star = o.elem.M
    v_planet = radvel(o)
    v_star = -(M_planet/(M_star + M_planet))*v_planet 
    return v_star
end
export radvel

"""
    acceleration(elem, t)

Get the instantaneous acceleration [mas/year^2] of
the *secondary* at the time `t` [days].

    acceleration(o)

Get the instantaneous acceleration [mas/year^2] of
the *secondary* from an instance of `AbstractOrbitSolution`.
"""
function acceleration(o::AbstractOrbitSolution)
    acc_planet = SVector(accra(o), accdec(o))
    return acc_planet
end

function accra(o::AbstractOrbitSolution)
    ẍcart = -o.elem.A*(1 + o.ecosν)^2 * (o.elem.cosi_cosΩ*o.sinν_ω + o.elem.sinΩ*o.cosν_ω) # [AU/year^2]
    cart2angle = rad2as*oftype(ẍcart, 1e3)/o.elem.dist
    ẍang = ẍcart*cart2angle # [mas/year^2] 
    return ẍang
end
function accdec(o::AbstractOrbitSolution)
    ÿcart = o.elem.A*(1 + o.ecosν)^2 * (o.elem.cosi_sinΩ*o.sinν_ω - o.elem.cosΩ*o.cosν_ω) # [AU/year^2]
    cart2angle = rad2as*oftype(ÿcart, 1e3)/o.elem.dist
    ÿang = ÿcart*cart2angle # [mas/year^2] 
    return ÿang
end
export accra, accdec

"""
    acceleration(elem, t, M_planet)

Get the instantaneous acceleration [mas/year^2] of 
the *primary* in at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.

acceleration(o)

Get the instantaneous acceleration [mas/year^2] of
the *primary* from an instance of `AbstractOrbitSolution`. The units of
`M_planet` and `elem.M` must match.
"""
function acceleration(o::AbstractOrbitSolution, M_planet)
    M_star = o.elem.M
    acc_planet = acceleration(o)
    acc_star = -(M_planet/(M_star + M_planet))*acc_planet
    return acc_star
end
export acceleration

# Define fallbacks for all accessor functions.
# If the user calls f(elems, t, args...) we compute the
# AbstractOrbitSolution for them.
fun_list = (
    :raoff,
    :decoff,
    :posangle,
    :projectedseparation,
    :propmotionanom,
    :radvel,
    :pmra,
    :pmdec,
    :accra,
    :accdec,
    :acceleration,
)
for fun in fun_list
    @eval function ($fun)(elem::AbstractOrbit, t::Real, args...)
        return ($fun)(orbitsolve(elem, t), args...)
    end
end

# ----------------------------------------------------------------------------------------------------------------------
# Kepler Equation Solver
# ----------------------------------------------------------------------------------------------------------------------

# The following function is taken directly from AstroLib.jl
# We  remove one invariant check we handle elsewhere and also
# force inlining for about a 5% speedup.
# We also supply analytic gradients for use in autodiff packages.
@inline function kepler_solver(_M::Real, e::Real)
    # We already handle this invariant
    # @assert 0 <= e <= 1 "eccentricity must be in the range [0, 1]"
    # M must be in the range [-pi, pi], see Markley (1995), page 2.
    M = rem2pi(_M, RoundNearest)
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



# ----------------------------------------------------------------------------------------------------------------------
# Addional & Optional Features
# ----------------------------------------------------------------------------------------------------------------------

include("diff-rules.jl")
include("transformation.jl")
include("time.jl")
include("plots-recipes.jl")

using Requires
function __init__()
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("makie-recipes.jl")

    # Small patch to allow symbolic tracing through the kepler solver.
    # Mean anomaly must still be in the range [0,2π] for the solution
    # to be valid.
    @require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" include("symbolics.jl")
end

end # module

# ----------------------------------------------------------------------------------------------------------------------